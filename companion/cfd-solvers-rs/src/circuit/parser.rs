//! SPICE netlist reader: engineering numbers, `{}` expressions with `.PARAM`, `.FUNC` macros,
//! `.SUBCKT` flattening with pins, default/instance parameters and local models, `.INCLUDE` /
//! `.LIB` file sections, element cards R C L K V I E F G H B D M Q J S (PULSE/SIN/PWL sources,
//! POLY(1) and TABLE dependent sources, behavioral `V()`/`TIME` expressions), `.MODEL` cards
//! and `.IC` / `.NODESET`. Ported from upstream `cfd/circuit/spice` (MIT, cfd_solvers
//! contributors).
#![allow(clippy::needless_range_loop)]
use super::spice::{
    upper, BehavioralEvaluator, BjtModel, Circuit, Complex, DiodeModel, JfetModel, MosLevel1Model,
    PulseWaveform, PwlWaveform, SineWaveform, SourceWaveform, SwitchModel,
};
use crate::workflow::text::{c_isspace, getline_split};
use crate::{require, Result, SolverError};
use std::{
    collections::{HashMap, HashSet},
    f64::consts::PI,
    path::{Path, PathBuf},
    sync::Arc,
};

fn invalid(msg: impl Into<String>) -> SolverError {
    SolverError::InvalidInput(msg.into())
}

fn is_space(c: u8) -> bool {
    c_isspace(c)
}

pub(crate) fn trim(s: &str) -> &str {
    s.trim_matches(|c: char| c.is_ascii() && is_space(c as u8))
}

pub(crate) fn tokens(line: &str) -> Vec<String> {
    line.split(|c: char| c.is_ascii() && is_space(c as u8))
        .filter(|t| !t.is_empty())
        .map(str::to_string)
        .collect()
}

fn ident_start(c: u8) -> bool {
    c.is_ascii_alphabetic() || c == b'_'
}

fn ident_char(c: u8) -> bool {
    c.is_ascii_alphanumeric() || c == b'_'
}

fn numeric_text(v: f64) -> String {
    crate::lbm::format_general(v, 17)
}

/// SPICE number with engineering suffix (T G MEG K M U N P F; MEG checked first).
pub fn parse_spice_number(text: &str) -> Result<f64> {
    let value = upper(trim(text));
    require(!value.is_empty(), "empty SPICE number")?;
    let (number, consumed) = stod_with_position(&value)?;
    let suffix = &value[consumed..];
    let scale = if suffix.starts_with("MEG") {
        1e6
    } else {
        match suffix.as_bytes().first() {
            Some(b'T') => 1e12,
            Some(b'G') => 1e9,
            Some(b'K') => 1e3,
            Some(b'M') => 1e-3,
            Some(b'U') => 1e-6,
            Some(b'N') => 1e-9,
            Some(b'P') => 1e-12,
            Some(b'F') => 1e-15,
            _ => 1.0,
        }
    };
    Ok(number * scale)
}

/// `std::stod` with the consumed length (longest numeric prefix).
fn stod_with_position(text: &str) -> Result<(f64, usize)> {
    for end in (1..=text.len()).rev() {
        if let Some(prefix) = text.get(..end) {
            if prefix.starts_with(|c: char| c.is_ascii_whitespace()) {
                break;
            }
            if let Ok(v) = prefix.parse::<f64>() {
                return Ok((v, end));
            }
        }
    }
    Err(invalid(format!("stod: no conversion for '{text}'")))
}

struct ExpressionParser<'a> {
    text: Vec<u8>,
    parameters: &'a HashMap<String, f64>,
    pos: usize,
}

impl<'a> ExpressionParser<'a> {
    fn new(text: &str, parameters: &'a HashMap<String, f64>) -> Self {
        let mut t = trim(text);
        if t.len() >= 2 && t.starts_with('{') && t.ends_with('}') {
            t = &t[1..t.len() - 1];
        }
        Self {
            text: t.as_bytes().to_vec(),
            parameters,
            pos: 0,
        }
    }
    fn parse(&mut self) -> Result<f64> {
        self.pos = 0;
        let v = self.expression()?;
        self.skip();
        require(
            self.pos == self.text.len(),
            "trailing SPICE expression text",
        )?;
        require(v.is_finite(), "non-finite SPICE expression")?;
        Ok(v)
    }
    fn skip(&mut self) {
        while self.pos < self.text.len() && is_space(self.text[self.pos]) {
            self.pos += 1;
        }
    }
    fn take(&mut self, c: u8) -> bool {
        self.skip();
        if self.pos < self.text.len() && self.text[self.pos] == c {
            self.pos += 1;
            true
        } else {
            false
        }
    }
    fn expression(&mut self) -> Result<f64> {
        let mut v = self.term()?;
        loop {
            if self.take(b'+') {
                v += self.term()?;
            } else if self.take(b'-') {
                v -= self.term()?;
            } else {
                return Ok(v);
            }
        }
    }
    fn term(&mut self) -> Result<f64> {
        let mut v = self.power()?;
        loop {
            if self.take(b'*') {
                v *= self.power()?;
            } else if self.take(b'/') {
                let d = self.power()?;
                require(d != 0.0, "division by zero in SPICE expression")?;
                v /= d;
            } else {
                return Ok(v);
            }
        }
    }
    fn power(&mut self) -> Result<f64> {
        let v = self.unary()?;
        if self.take(b'^') {
            return Ok(v.powf(self.power()?));
        }
        Ok(v)
    }
    fn unary(&mut self) -> Result<f64> {
        if self.take(b'+') {
            return self.unary();
        }
        if self.take(b'-') {
            return Ok(-self.unary()?);
        }
        self.primary()
    }
    fn primary(&mut self) -> Result<f64> {
        self.skip();
        if self.take(b'(') {
            let v = self.expression()?;
            require(self.take(b')'), "missing ) in SPICE expression")?;
            return Ok(v);
        }
        require(
            self.pos < self.text.len(),
            "unexpected end of SPICE expression",
        )?;
        if ident_start(self.text[self.pos]) {
            let begin = self.pos;
            self.pos += 1;
            while self.pos < self.text.len() && ident_char(self.text[self.pos]) {
                self.pos += 1;
            }
            let key = upper(&String::from_utf8_lossy(&self.text[begin..self.pos]));
            return self
                .parameters
                .get(&key)
                .copied()
                .ok_or_else(|| invalid(format!("unknown SPICE parameter: {key}")));
        }
        let begin = self.pos;
        let mut exponent = false;
        while self.pos < self.text.len() {
            let c = self.text[self.pos];
            if c.is_ascii_digit() || c == b'.' {
                self.pos += 1;
                continue;
            }
            if (c == b'e' || c == b'E') && !exponent {
                exponent = true;
                self.pos += 1;
                if self.pos < self.text.len()
                    && (self.text[self.pos] == b'+' || self.text[self.pos] == b'-')
                {
                    self.pos += 1;
                }
                continue;
            }
            break;
        }
        require(begin != self.pos, "expected SPICE expression value")?;
        while self.pos < self.text.len() && self.text[self.pos].is_ascii_alphabetic() {
            self.pos += 1;
        }
        parse_spice_number(&String::from_utf8_lossy(&self.text[begin..self.pos]))
    }
}

/// Evaluate a SPICE expression (`+ - * / ^`, parentheses, unary signs, engineering numbers and
/// case-insensitive parameter names; optional surrounding `{}`).
pub fn evaluate_spice_expression(
    expression: &str,
    parameters: &HashMap<String, f64>,
) -> Result<f64> {
    let normalized: HashMap<String, f64> = parameters.iter().map(|(k, v)| (upper(k), *v)).collect();
    ExpressionParser::new(expression, &normalized).parse()
}

fn eval_upper(expression: &str, parameters: &HashMap<String, f64>) -> Result<f64> {
    ExpressionParser::new(expression, parameters).parse()
}

struct UserFunction {
    name: String,
    arguments: Vec<String>,
    expression: String,
}

fn replace_identifier_ci(text: &str, identifier: &str, replacement: &str) -> String {
    let key = upper(identifier);
    let b = text.as_bytes();
    let mut out: Vec<u8> = Vec::with_capacity(text.len() + replacement.len());
    let mut i = 0;
    while i < b.len() {
        if ident_start(b[i]) {
            let begin = i;
            i += 1;
            while i < b.len() && ident_char(b[i]) {
                i += 1;
            }
            let token = &text[begin..i];
            if upper(token) == key {
                out.extend_from_slice(replacement.as_bytes());
            } else {
                out.extend_from_slice(token.as_bytes());
            }
        } else {
            out.push(b[i]);
            i += 1;
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn split_function_arguments(text: &str) -> Vec<String> {
    let b = text.as_bytes();
    let mut out = Vec::new();
    let (mut begin, mut depth) = (0usize, 0i32);
    for i in 0..=b.len() {
        let c = if i < b.len() { b[i] } else { b',' };
        if c == b'(' || c == b'{' {
            depth += 1;
        } else if c == b')' || c == b'}' {
            depth -= 1;
        }
        if c == b',' && depth == 0 {
            out.push(trim(&text[begin..i]).to_string());
            begin = i + 1;
        }
    }
    if out.len() == 1 && out[0].is_empty() {
        out.clear();
    }
    out
}

fn expand_user_functions(
    text: &str,
    functions: &HashMap<String, UserFunction>,
    depth: usize,
) -> Result<String> {
    require(depth <= 32, "SPICE .FUNC expansion depth exceeded")?;
    let mut text = text.to_string();
    let mut i = 0;
    while i < text.len() {
        let b = text.as_bytes();
        if !ident_start(b[i]) {
            i += 1;
            continue;
        }
        let begin = i;
        i += 1;
        while i < b.len() && ident_char(b[i]) {
            i += 1;
        }
        let name = upper(&text[begin..i]);
        let mut open = i;
        while open < b.len() && is_space(b[open]) {
            open += 1;
        }
        let Some(function) = functions.get(&name) else {
            continue;
        };
        if open >= b.len() || b[open] != b'(' {
            continue;
        }
        let mut nesting = 1;
        let mut close = open + 1;
        while close < b.len() && nesting > 0 {
            if b[close] == b'(' {
                nesting += 1;
            } else if b[close] == b')' {
                nesting -= 1;
            }
            close += 1;
        }
        require(nesting == 0, "unterminated SPICE .FUNC call")?;
        let close_index = close - 1;
        let actual = split_function_arguments(&text[open + 1..close_index]);
        require(
            actual.len() == function.arguments.len(),
            format!("SPICE .FUNC argument count mismatch: {}", function.name),
        )?;
        let mut expression = function.expression.clone();
        for (a, arg) in actual.iter().enumerate() {
            let expanded = expand_user_functions(arg, functions, depth + 1)?;
            expression = replace_identifier_ci(
                &expression,
                &function.arguments[a],
                &format!("({expanded})"),
            );
        }
        expression = expand_user_functions(&expression, functions, depth + 1)?;
        text.replace_range(begin..=close_index, &format!("({expression})"));
        i = begin;
    }
    Ok(text)
}

fn preprocess_user_functions(lines: &[String]) -> Result<Vec<String>> {
    let mut functions: HashMap<String, UserFunction> = HashMap::new();
    let mut body = Vec::with_capacity(lines.len());
    for raw in lines {
        let line = trim(raw);
        let t = tokens(line);
        if t.is_empty() || upper(&t[0]) != ".FUNC" {
            body.push(raw.clone());
            continue;
        }
        let name_start = line
            .bytes()
            .enumerate()
            .skip(5)
            .find(|&(_, c)| c != b' ' && c != b'\t')
            .map(|(i, _)| i)
            .ok_or_else(|| invalid("SPICE .FUNC missing definition"))?;
        let open = line[name_start..].find('(').map(|o| o + name_start);
        let close = open.and_then(|o| line[o + 1..].find(')').map(|c| c + o + 1));
        let (Some(open), Some(close)) = (open, close) else {
            return Err(invalid("malformed SPICE .FUNC signature"));
        };
        let name = upper(trim(&line[name_start..open]));
        let arguments = split_function_arguments(&line[open + 1..close])
            .iter()
            .map(|a| upper(a))
            .collect();
        let mut expression = trim(&line[close + 1..]).to_string();
        if expression.len() >= 2 && expression.starts_with('{') && expression.ends_with('}') {
            expression = expression[1..expression.len() - 1].to_string();
        }
        require(
            !name.is_empty() && !expression.is_empty(),
            "malformed SPICE .FUNC definition",
        )?;
        functions.insert(
            name.clone(),
            UserFunction {
                name,
                arguments,
                expression,
            },
        );
    }
    if functions.is_empty() {
        return Ok(body);
    }
    body.iter()
        .map(|l| expand_user_functions(l, &functions, 0))
        .collect()
}

fn unquote_path(token: &str) -> String {
    let t = trim(token);
    if t.len() >= 2
        && ((t.starts_with('"') && t.ends_with('"')) || (t.starts_with('\'') && t.ends_with('\'')))
    {
        return t[1..t.len() - 1].to_string();
    }
    t.to_string()
}

fn load_spice_file_recursive(
    input: &Path,
    section: &str,
    active: &mut HashSet<String>,
) -> Result<Vec<String>> {
    let path = std::fs::canonicalize(input).unwrap_or_else(|_| input.to_path_buf());
    let key = path.to_string_lossy().into_owned();
    require(
        active.insert(key.clone()),
        format!("recursive SPICE include: {key}"),
    )?;
    let bytes = std::fs::read(&path)
        .map_err(|e| SolverError::Io(format!("cannot open SPICE file: {key}: {e}")))?;
    let text = String::from_utf8_lossy(&bytes).into_owned();
    let base = path.parent().map(Path::to_path_buf).unwrap_or_default();
    let mut selected = Vec::new();
    let mut section_active = section.is_empty();
    let mut found_section = section.is_empty();
    let mut library_depth = 0usize;
    for raw in getline_split(&text, b'\n') {
        let t = tokens(trim(raw));
        if t.is_empty() {
            if section_active {
                selected.push(raw.to_string());
            }
            continue;
        }
        let directive = upper(&t[0]);
        if !section.is_empty() {
            if directive == ".LIB" && t.len() == 2 && upper(&t[1]) == upper(section) {
                section_active = true;
                found_section = true;
                continue;
            }
            if section_active && directive == ".ENDL" {
                section_active = false;
                continue;
            }
            if !section_active {
                continue;
            }
        } else {
            if directive == ".LIB" && t.len() == 2 {
                library_depth += 1;
                continue;
            }
            if directive == ".ENDL" && library_depth > 0 {
                library_depth -= 1;
                continue;
            }
            if library_depth > 0 {
                continue;
            }
        }
        if directive == ".INCLUDE" {
            require(t.len() >= 2, "SPICE .INCLUDE requires a path")?;
            selected.extend(load_spice_file_recursive(
                &base.join(unquote_path(&t[1])),
                "",
                active,
            )?);
            continue;
        }
        if directive == ".LIB" && t.len() >= 3 {
            selected.extend(load_spice_file_recursive(
                &base.join(unquote_path(&t[1])),
                &t[2],
                active,
            )?);
            continue;
        }
        selected.push(raw.to_string());
    }
    active.remove(&key);
    require(
        found_section,
        format!("SPICE .LIB section not found: {section}"),
    )?;
    Ok(selected)
}

struct SubcircuitDefinition {
    name: String,
    pins: Vec<String>,
    default_parameters: Vec<(String, String)>,
    body: Vec<String>,
}

fn parse_parameter_assignments(
    assignments: &[String],
    inherited: &HashMap<String, f64>,
) -> Result<HashMap<String, f64>> {
    let mut parameters = inherited.clone();
    for token in assignments {
        if upper(token) == "PARAMS:" {
            continue;
        }
        let Some(eq) = token.find('=') else {
            continue;
        };
        let key = upper(&token[..eq]);
        let v = eval_upper(&token[eq + 1..], &parameters)?;
        parameters.insert(key, v);
    }
    Ok(parameters)
}

fn map_subcircuit_node(token: &str, node_map: &HashMap<String, String>, prefix: &str) -> String {
    if token == "0" {
        return token.to_string();
    }
    match node_map.get(&upper(token)) {
        Some(m) => m.clone(),
        None => format!("{prefix}{token}"),
    }
}

fn substitute_parameter_token(token: &str, parameters: &HashMap<String, f64>) -> Result<String> {
    if token.is_empty() {
        return Ok(String::new());
    }
    if let Some(v) = parameters.get(&upper(token)) {
        return Ok(numeric_text(*v));
    }
    let mut token = token.to_string();
    let mut begin = 0;
    while let Some(off) = token.get(begin..).and_then(|s| s.find('{')) {
        let start = begin + off;
        let end = token[start + 1..]
            .find('}')
            .map(|e| e + start + 1)
            .ok_or_else(|| invalid("unclosed parameter expression in subcircuit"))?;
        let v = eval_upper(&token[start..=end], parameters)?;
        token.replace_range(start..=end, &numeric_text(v));
        begin = start + 1;
    }
    Ok(token)
}

/// Replace every `V(a[,b])` (case-insensitive V) by applying `f` to the argument list.
fn rewrite_voltage_references(
    text: &str,
    mut f: impl FnMut(&str, &[String], &mut String) -> Result<()>,
) -> Result<String> {
    let b = text.as_bytes();
    let mut out: Vec<u8> = Vec::with_capacity(text.len() + 16);
    let mut i = 0;
    while i < b.len() {
        if (b[i] == b'V' || b[i] == b'v') && i + 1 < b.len() && b[i + 1] == b'(' {
            let close = text[i + 2..]
                .find(')')
                .map(|c| c + i + 2)
                .ok_or_else(|| invalid("unclosed V() reference in behavioral source"))?;
            let args = split_function_arguments(&text[i + 2..close]);
            let mut piece = String::new();
            f(&text[i..=close], &args, &mut piece)?;
            out.extend_from_slice(piece.as_bytes());
            i = close + 1;
            continue;
        }
        out.push(b[i]);
        i += 1;
    }
    Ok(String::from_utf8_lossy(&out).into_owned())
}

fn map_behavior_voltage_references(
    text: &str,
    node_map: &HashMap<String, String>,
    prefix: &str,
) -> Result<String> {
    rewrite_voltage_references(text, |_, args, out| {
        require(
            !args.is_empty() && args.len() <= 2,
            "behavioral V() requires one or two node names",
        )?;
        out.push_str("V(");
        out.push_str(&map_subcircuit_node(&args[0], node_map, prefix));
        if args.len() == 2 {
            out.push(',');
            out.push_str(&map_subcircuit_node(&args[1], node_map, prefix));
        }
        out.push(')');
        Ok(())
    })
}

fn substitute_behavior_parameters(text: &str, parameters: &HashMap<String, f64>) -> Result<String> {
    let b = text.as_bytes();
    let mut out: Vec<u8> = Vec::with_capacity(text.len() + 16);
    let mut i = 0;
    while i < b.len() {
        if (b[i] == b'V' || b[i] == b'v') && i + 1 < b.len() && b[i + 1] == b'(' {
            let close = text[i + 2..]
                .find(')')
                .map(|c| c + i + 2)
                .ok_or_else(|| invalid("unclosed V() reference in behavioral source"))?;
            out.extend_from_slice(&b[i..=close]);
            i = close + 1;
            continue;
        }
        if ident_start(b[i]) {
            let begin = i;
            i += 1;
            while i < b.len() && ident_char(b[i]) {
                i += 1;
            }
            let token = &text[begin..i];
            match parameters.get(&upper(token)) {
                Some(v) => out.extend_from_slice(numeric_text(*v).as_bytes()),
                None => out.extend_from_slice(token.as_bytes()),
            }
            continue;
        }
        out.push(b[i]);
        i += 1;
    }
    Ok(String::from_utf8_lossy(&out).into_owned())
}

fn first_upper(t: &str) -> u8 {
    t.as_bytes().first().map_or(0, |c| c.to_ascii_uppercase())
}

fn flatten_device_line(
    raw: &str,
    node_map: &HashMap<String, String>,
    parameters: &HashMap<String, f64>,
    model_map: &HashMap<String, String>,
    prefix: &str,
) -> Result<String> {
    let mut t = tokens(raw);
    if t.is_empty() {
        return Ok(String::new());
    }
    let kind = first_upper(&t[0]);
    let prefixed = |name: &str| -> String {
        if name.is_empty() || prefix.is_empty() {
            name.to_string()
        } else {
            let split = name.chars().next().map_or(0, char::len_utf8);
            format!("{}{prefix}{}", &name[..split], &name[split..])
        }
    };
    t[0] = prefixed(&t[0]);
    let map_nodes = |t: &mut Vec<String>, count: usize| {
        for i in 1..=count.min(t.len().saturating_sub(1)) {
            t[i] = map_subcircuit_node(&t[i], node_map, prefix);
        }
    };
    let dependent_mode = if matches!(kind, b'G' | b'E') && t.len() > 3 {
        upper(&t[3])
    } else {
        String::new()
    };
    let poly = dependent_mode.starts_with("POLY(");
    let table = dependent_mode == "TABLE";
    match kind {
        b'R' | b'C' | b'L' | b'D' | b'V' | b'I' | b'B' => map_nodes(&mut t, 2),
        b'G' | b'E' if poly || table => {
            map_nodes(&mut t, 2);
            if poly && t.len() >= 6 {
                t[4] = map_subcircuit_node(&t[4], node_map, prefix);
                t[5] = map_subcircuit_node(&t[5], node_map, prefix);
            }
        }
        b'G' | b'E' | b'S' => map_nodes(&mut t, 4),
        b'F' | b'H' => map_nodes(&mut t, 2),
        b'M' => map_nodes(&mut t, 4),
        b'Q' | b'J' => map_nodes(&mut t, 3),
        _ => {}
    }
    if kind == b'K' && t.len() >= 4 {
        t[1] = prefixed(&t[1]);
        t[2] = prefixed(&t[2]);
    }
    if matches!(kind, b'F' | b'H') && t.len() >= 4 {
        t[3] = prefixed(&t[3]);
    }
    let model_index = match kind {
        b'D' => Some(3),
        b'M' | b'S' => Some(5),
        b'Q' | b'J' => Some(4),
        _ => None,
    };
    if let Some(mi) = model_index.filter(|&mi| mi < t.len()) {
        if let Some(m) = model_map.get(&upper(&t[mi])) {
            t[mi] = m.clone();
        }
    }
    let join_from = |t: &[String], start: usize| t[start..].join(" ");
    if matches!(kind, b'R' | b'C' | b'L') && t.len() >= 4 {
        t[3] = substitute_parameter_token(&t[3], parameters)?;
    } else if matches!(kind, b'G' | b'E') && poly {
        for i in 6..t.len() {
            t[i] = substitute_parameter_token(&t[i], parameters)?;
        }
    } else if matches!(kind, b'G' | b'E') && table {
        let spec = map_behavior_voltage_references(&join_from(&t, 4), node_map, prefix)?;
        let spec = substitute_behavior_parameters(&spec, parameters)?;
        t.truncate(5);
        t[4] = spec;
    } else if matches!(kind, b'G' | b'E') && t.len() >= 6 {
        t[5] = substitute_parameter_token(&t[5], parameters)?;
    } else if matches!(kind, b'F' | b'H') && t.len() >= 5 {
        t[4] = substitute_parameter_token(&t[4], parameters)?;
    } else if kind == b'K' && t.len() >= 4 {
        t[3] = substitute_parameter_token(&t[3], parameters)?;
    } else if matches!(kind, b'V' | b'I') {
        for i in 3..t.len() {
            t[i] = substitute_parameter_token(&t[i], parameters)?;
        }
    } else if kind == b'B' && t.len() >= 4 {
        let spec = map_behavior_voltage_references(&join_from(&t, 3), node_map, prefix)?;
        let spec = substitute_behavior_parameters(&spec, parameters)?;
        t.truncate(4);
        t[3] = spec;
    }
    Ok(t.join(" "))
}

struct Flattener<'a> {
    definitions: &'a [SubcircuitDefinition],
    output: Vec<String>,
    stack: Vec<String>,
}

impl Flattener<'_> {
    fn expand(
        &mut self,
        instance: &[String],
        inherited_parameters: &HashMap<String, f64>,
        outer_nodes: &HashMap<String, String>,
        inherited_models: &HashMap<String, String>,
        outer_prefix: &str,
    ) -> Result<()> {
        let mut found = None;
        for d in self.definitions {
            let index = 1 + d.pins.len();
            if index < instance.len() && upper(&instance[index]) == d.name {
                found = Some((d, index));
                break;
            }
        }
        let Some((definition, name_index)) = found else {
            return Err(invalid(format!(
                "unknown or malformed SPICE subcircuit instance: {}",
                instance.first().map_or("", String::as_str)
            )));
        };
        require(
            !self.stack.contains(&definition.name),
            "recursive SPICE subcircuit hierarchy",
        )?;
        let mut parameters = inherited_parameters.clone();
        for (key, expression) in &definition.default_parameters {
            let v = eval_upper(expression, &parameters)?;
            parameters.insert(key.clone(), v);
        }
        parameters = parse_parameter_assignments(&instance[name_index + 1..], &parameters)?;
        let mut node_map = HashMap::new();
        for (i, pin) in definition.pins.iter().enumerate() {
            node_map.insert(
                pin.clone(),
                map_subcircuit_node(&instance[i + 1], outer_nodes, outer_prefix),
            );
        }
        let prefix = format!("{outer_prefix}{}:", instance[0]);
        let mut models = inherited_models.clone();
        for line in &definition.body {
            let mt = tokens(line);
            if mt.len() >= 3 && upper(&mt[0]) == ".MODEL" {
                models.insert(upper(&mt[1]), format!("{prefix}{}", mt[1]));
            }
        }
        self.stack.push(definition.name.clone());
        for line in &definition.body {
            let mut t = tokens(line);
            if t.is_empty() {
                continue;
            }
            let head = upper(&t[0]);
            if head == ".PARAM" {
                parameters = parse_parameter_assignments(&t[1..], &parameters)?;
                continue;
            }
            if head == ".MODEL" {
                require(t.len() >= 3, "malformed local SPICE .MODEL")?;
                t[1] = models
                    .get(&upper(&t[1]))
                    .cloned()
                    .ok_or_else(|| invalid("unknown local SPICE model"))?;
                for k in 2..t.len() {
                    t[k] = substitute_parameter_token(&t[k], &parameters)?;
                }
                self.output.push(t.join(" "));
                continue;
            }
            if first_upper(&t[0]) == b'X' {
                self.expand(&t, &parameters, &node_map, &models, &prefix)?;
            } else {
                let flat = flatten_device_line(line, &node_map, &parameters, &models, &prefix)?;
                self.output.push(flat);
            }
        }
        self.stack.pop();
        Ok(())
    }
}

fn flatten_subcircuits(lines: &[String]) -> Result<Vec<String>> {
    let mut definitions: Vec<SubcircuitDefinition> = Vec::new();
    let mut top_level = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        let line = trim(&lines[i]).to_string();
        if line.is_empty() || line.starts_with('*') {
            i += 1;
            continue;
        }
        let t = tokens(&line);
        if !t.is_empty() && upper(&t[0]) == ".SUBCKT" {
            require(t.len() >= 3, "SPICE .SUBCKT requires name and pins")?;
            let name = upper(&t[1]);
            let parameter_start = (2..t.len())
                .find(|&k| upper(&t[k]) == "PARAMS:" || t[k].contains('='))
                .unwrap_or(t.len());
            let pins = t[2..parameter_start].iter().map(|p| upper(p)).collect();
            let mut default_parameters: Vec<(String, String)> = Vec::new();
            for tok in &t[parameter_start..] {
                if upper(tok) == "PARAMS:" {
                    continue;
                }
                if let Some(eq) = tok.find('=') {
                    let key = upper(&tok[..eq]);
                    let value = tok[eq + 1..].to_string();
                    match default_parameters.iter_mut().find(|(k, _)| *k == key) {
                        Some(slot) => slot.1 = value,
                        None => default_parameters.push((key, value)),
                    }
                }
            }
            let mut body = Vec::new();
            let mut ended = false;
            i += 1;
            while i < lines.len() {
                let b = trim(&lines[i]).to_string();
                let bt = tokens(&b);
                if !bt.is_empty() && upper(&bt[0]) == ".ENDS" {
                    ended = true;
                    break;
                }
                body.push(b);
                i += 1;
            }
            require(ended, format!("unterminated SPICE .SUBCKT {name}"))?;
            let def = SubcircuitDefinition {
                name,
                pins,
                default_parameters,
                body,
            };
            match definitions.iter_mut().find(|d| d.name == def.name) {
                Some(slot) => *slot = def,
                None => definitions.push(def),
            }
        } else {
            top_level.push(line);
        }
        i += 1;
    }
    if definitions.is_empty() {
        return Ok(top_level);
    }
    let mut global = HashMap::new();
    for line in &top_level {
        let t = tokens(line);
        if !t.is_empty() && upper(&t[0]) == ".PARAM" {
            global = parse_parameter_assignments(&t[1..], &global)?;
        }
    }
    let mut flattener = Flattener {
        definitions: &definitions,
        output: Vec::new(),
        stack: Vec::new(),
    };
    let none = HashMap::new();
    for line in &top_level {
        let t = tokens(line);
        if t.is_empty() {
            continue;
        }
        if first_upper(&t[0]) == b'X' {
            flattener.expand(&t, &global, &none, &none, "")?;
        } else {
            flattener.output.push(line.clone());
        }
    }
    Ok(flattener.output)
}

fn model_params(line: &str) -> HashMap<String, String> {
    let mapped: String = line
        .chars()
        .map(|c| if matches!(c, '(' | ')' | ',') { ' ' } else { c })
        .collect();
    let mut out = HashMap::new();
    for t in tokens(&mapped) {
        if let Some(eq) = t.find('=') {
            out.insert(upper(&t[..eq]), t[eq + 1..].to_string());
        }
    }
    out
}

struct VoltageReference {
    positive: usize,
    negative: usize,
    parameter: String,
}

fn behavioral_evaluator(
    circuit: &mut Circuit,
    expression: &str,
    parameters: &HashMap<String, f64>,
) -> Result<BehavioralEvaluator> {
    let expression = trim(expression);
    require(!expression.is_empty(), "empty behavioral-source expression")?;
    let mut references: Vec<VoltageReference> = Vec::new();
    let transformed = rewrite_voltage_references(expression, |_, args, out| {
        require(
            !args.is_empty() && args.len() <= 2,
            "behavioral V() requires one or two node names",
        )?;
        let parameter = format!("__BV{}", references.len());
        references.push(VoltageReference {
            positive: circuit.node(&args[0]),
            negative: if args.len() == 2 {
                circuit.node(&args[1])
            } else {
                0
            },
            parameter: parameter.clone(),
        });
        out.push_str(&parameter);
        Ok(())
    })?;
    let captured = parameters.clone();
    Ok(Arc::new(move |nodes: &[f64], time: f64| {
        let mut local = captured.clone();
        local.insert("TIME".into(), time);
        for r in &references {
            let (Some(p), Some(n)) = (nodes.get(r.positive), nodes.get(r.negative)) else {
                return Err(invalid("behavioral-source node reference out of range"));
            };
            local.insert(r.parameter.clone(), p - n);
        }
        eval_upper(&transformed, &local)
    }))
}

fn source_waveform(
    specification: &str,
    value: &dyn Fn(&str) -> Result<f64>,
) -> Result<(f64, SourceWaveform)> {
    let mapped: String = specification
        .chars()
        .map(|c| if matches!(c, '(' | ')' | ',') { ' ' } else { c })
        .collect();
    let args = tokens(&mapped);
    let kind = upper(
        args.first()
            .ok_or_else(|| invalid("empty transient source specification"))?,
    );
    match kind.as_str() {
        "PULSE" => {
            require(args.len() >= 8, "PULSE requires V1 V2 TD TR TF PW PER")?;
            let p = PulseWaveform {
                initial: value(&args[1])?,
                pulsed: value(&args[2])?,
                delay: value(&args[3])?,
                rise: value(&args[4])?,
                fall: value(&args[5])?,
                width: value(&args[6])?,
                period: value(&args[7])?,
            };
            Ok((p.initial, p.waveform()))
        }
        "SIN" => {
            require(args.len() >= 4, "SIN requires VO VA FREQ")?;
            let mut s = SineWaveform {
                offset: value(&args[1])?,
                amplitude: value(&args[2])?,
                frequency_hz: value(&args[3])?,
                ..Default::default()
            };
            let opt = |i: usize| args.get(i).filter(|a| upper(a) != "AC");
            if let Some(a) = opt(4) {
                s.delay = value(a)?;
            }
            if let Some(a) = opt(5) {
                s.damping_per_s = value(a)?;
            }
            if let Some(a) = opt(6) {
                s.phase_deg = value(a)?;
            }
            Ok((s.offset, s.waveform()))
        }
        "PWL" => {
            let mut pwl = PwlWaveform::default();
            let mut i = 1;
            while i + 1 < args.len() && upper(&args[i]) != "AC" {
                pwl.points.push((value(&args[i])?, value(&args[i + 1])?));
                i += 2;
            }
            let dc = pwl
                .points
                .first()
                .map(|p| p.1)
                .ok_or_else(|| invalid("PWL requires time/value pairs"))?;
            Ok((dc, pwl.waveform()))
        }
        _ => Err(invalid("unsupported transient source waveform")),
    }
}

fn apply_voltage_directive(
    circuit: &mut Circuit,
    raw: &str,
    initial: bool,
    value: &dyn Fn(&str) -> Result<f64>,
) -> Result<()> {
    let first_space = raw
        .find([' ', '\t'])
        .ok_or_else(|| invalid("empty SPICE voltage directive"))?;
    let body = &raw.as_bytes()[first_space + 1..];
    let text = &raw[first_space + 1..];
    let mut pos = 0;
    let mut any = false;
    while pos < body.len() {
        while pos < body.len() && is_space(body[pos]) {
            pos += 1;
        }
        if pos >= body.len() {
            break;
        }
        require(
            body[pos].eq_ignore_ascii_case(&b'V') && pos + 1 < body.len() && body[pos + 1] == b'(',
            "SPICE .IC/.NODESET currently expects V(node)=value assignments",
        )?;
        let close = text[pos + 2..]
            .find(')')
            .map(|c| c + pos + 2)
            .ok_or_else(|| invalid("malformed SPICE voltage directive node"))?;
        let node_name = trim(&text[pos + 2..close]);
        require(
            !node_name.is_empty() && !node_name.contains(','),
            "SPICE .IC/.NODESET differential assignments are not supported yet",
        )?;
        let mut eq = close + 1;
        while eq < body.len() && is_space(body[eq]) {
            eq += 1;
        }
        require(
            eq < body.len() && body[eq] == b'=',
            "SPICE voltage directive requires =",
        )?;
        let mut begin = eq + 1;
        while begin < body.len() && is_space(body[begin]) {
            begin += 1;
        }
        let mut end = begin;
        let mut depth = 0i32;
        while end < body.len() {
            let c = body[end];
            if c == b'{' || c == b'(' {
                depth += 1;
            } else if c == b'}' || c == b')' {
                depth -= 1;
            }
            if depth == 0 && is_space(c) {
                break;
            }
            end += 1;
        }
        require(begin != end, "SPICE voltage directive missing value")?;
        let voltage = value(&text[begin..end])?;
        let node = circuit.node(node_name);
        if initial {
            circuit.set_initial_voltage(node, voltage)?;
        } else {
            circuit.set_nodeset_voltage(node, voltage)?;
        }
        any = true;
        pos = end;
    }
    require(any, "empty SPICE voltage directive")
}

fn joined(t: &[String], start: usize) -> String {
    t[start..].join(" ")
}

impl Circuit {
    /// Parse a netlist held in memory (`.INCLUDE` / `.LIB file section` need
    /// [`Circuit::parse_spice_file`]).
    pub fn parse_spice(text: &str) -> Result<Circuit> {
        let raw_lines: Vec<String> = getline_split(text, b'\n')
            .into_iter()
            .map(str::to_string)
            .collect();
        for raw in &raw_lines {
            let t = tokens(trim(raw));
            if !t.is_empty()
                && (upper(&t[0]) == ".INCLUDE" || (upper(&t[0]) == ".LIB" && t.len() >= 3))
            {
                return Err(invalid("file directives require Circuit::parse_spice_file"));
            }
        }
        let functional = preprocess_user_functions(&raw_lines)?;
        let flattened = flatten_subcircuits(&functional)?;
        let (mut devices, mut mutuals, mut models, mut params, mut ics, mut nodesets) =
            (vec![], vec![], vec![], vec![], vec![], vec![]);
        for line in &flattened {
            let line = trim(line).to_string();
            if line.is_empty() || line.starts_with('*') {
                continue;
            }
            if line.starts_with('.') {
                let head = tokens(&line).first().map(|h| upper(h)).unwrap_or_default();
                match head.as_str() {
                    ".MODEL" => models.push(line),
                    ".PARAM" => params.push(line),
                    ".IC" => ics.push(line),
                    ".NODESET" => nodesets.push(line),
                    _ => {}
                }
                continue;
            }
            if first_upper(&line) == b'K' {
                mutuals.push(line);
            } else {
                devices.push(line);
            }
        }
        let mut parameters: HashMap<String, f64> = HashMap::new();
        for raw in &params {
            let t = tokens(&raw.replace(',', " "));
            for tok in t.iter().skip(1) {
                let Some(eq) = tok.find('=') else {
                    continue;
                };
                let v = eval_upper(&tok[eq + 1..], &parameters)?;
                parameters.insert(upper(&tok[..eq]), v);
            }
        }
        let value = |token: &str| eval_upper(token, &parameters);
        let mut circuit = Circuit::new();

        for raw in &models {
            let normalized: String = raw
                .chars()
                .map(|c| if matches!(c, '(' | ')' | ',') { ' ' } else { c })
                .collect();
            let t = tokens(&normalized);
            if t.len() < 3 {
                continue;
            }
            let (name, kind) = (t[1].as_str(), upper(&t[2]));
            let mp = model_params(raw);
            let mv = |key: &str, fallback: f64| -> Result<f64> {
                match mp.get(key) {
                    Some(v) => value(v),
                    None => Ok(fallback),
                }
            };
            match kind.as_str() {
                "D" => {
                    let d = DiodeModel::default();
                    circuit.add_diode_model(
                        name,
                        DiodeModel {
                            saturation_current: mv("IS", d.saturation_current)?,
                            emission_coefficient: mv("N", d.emission_coefficient)?,
                            ..d
                        },
                    )?;
                }
                "NMOS" | "PMOS" => {
                    let d = MosLevel1Model::default();
                    circuit.add_mos_model(
                        name,
                        MosLevel1Model {
                            pmos: kind == "PMOS",
                            threshold_voltage: mv("VTO", d.threshold_voltage)?,
                            transconductance_parameter: mv("KP", d.transconductance_parameter)?,
                            channel_length_modulation: mv("LAMBDA", d.channel_length_modulation)?,
                        },
                    )?;
                }
                "NPN" | "PNP" => {
                    let d = BjtModel::default();
                    circuit.add_bjt_model(
                        name,
                        BjtModel {
                            pnp: kind == "PNP",
                            saturation_current: mv("IS", d.saturation_current)?,
                            forward_beta: mv("BF", d.forward_beta)?,
                            reverse_beta: mv("BR", d.reverse_beta)?,
                            ..d
                        },
                    )?;
                }
                "NJF" | "PJF" => {
                    let d = JfetModel::default();
                    circuit.add_jfet_model(
                        name,
                        JfetModel {
                            pchannel: kind == "PJF",
                            pinch_off_voltage: mv("VTO", d.pinch_off_voltage)?,
                            beta: mv("BETA", d.beta)?,
                            channel_length_modulation: mv("LAMBDA", d.channel_length_modulation)?,
                        },
                    )?;
                }
                "SW" => {
                    let d = SwitchModel::default();
                    circuit.add_switch_model(
                        name,
                        SwitchModel {
                            on_resistance: mv("RON", d.on_resistance)?,
                            off_resistance: mv("ROFF", d.off_resistance)?,
                            threshold_voltage: mv("VT", d.threshold_voltage)?,
                            hysteresis_voltage: mv("VH", d.hysteresis_voltage)?,
                        },
                    )?;
                }
                _ => {}
            }
        }

        for raw in &devices {
            let t = tokens(raw);
            if t.is_empty() {
                continue;
            }
            let kind = first_upper(&t[0]);
            match kind {
                b'R' | b'C' | b'L' if t.len() >= 4 => {
                    let (p, n) = (circuit.node(&t[1]), circuit.node(&t[2]));
                    let v = value(&t[3])?;
                    match kind {
                        b'R' => circuit.add_resistor(t[0].clone(), p, n, v)?,
                        b'C' => circuit.add_capacitor(t[0].clone(), p, n, v)?,
                        _ => circuit.add_inductor(t[0].clone(), p, n, v)?,
                    }
                }
                b'V' | b'I' if t.len() >= 4 => {
                    let (p, n) = (circuit.node(&t[1]), circuit.node(&t[2]));
                    let dc;
                    let mut ac = Complex::default();
                    let mut transient = None;
                    let first = upper(&t[3]);
                    if first.starts_with("PULSE")
                        || first.starts_with("SIN")
                        || first.starts_with("PWL")
                    {
                        let (d, w) = source_waveform(&joined(&t, 3), &value)?;
                        dc = d;
                        transient = Some(w);
                    } else if upper(&t[3]) == "DC" && 4 < t.len() {
                        dc = value(&t[4])?;
                    } else {
                        dc = value(&t[3])?;
                    }
                    let mut i = 3;
                    while i < t.len() {
                        if upper(&t[i]) == "AC" && i + 1 < t.len() {
                            i += 1;
                            let magnitude = value(&t[i])?;
                            let mut phase = 0.0;
                            if i + 1 < t.len() {
                                i += 1;
                                phase = value(&t[i])?;
                            }
                            ac = Complex::from_polar(magnitude, phase * PI / 180.0);
                            break;
                        }
                        i += 1;
                    }
                    if kind == b'V' {
                        circuit.add_voltage_source(t[0].clone(), p, n, dc, ac, transient)?;
                    } else {
                        circuit.add_current_source(t[0].clone(), p, n, dc, ac, transient)?;
                    }
                }
                b'B' if t.len() >= 4 => {
                    let (p, n) = (circuit.node(&t[1]), circuit.node(&t[2]));
                    let spec = joined(&t, 3);
                    let eq = spec.find('=').ok_or_else(|| {
                        invalid("behavioral source requires V={expr} or I={expr}")
                    })?;
                    let output = upper(trim(&spec[..eq]));
                    let evaluator =
                        behavioral_evaluator(&mut circuit, &spec[eq + 1..], &parameters)?;
                    match output.as_str() {
                        "V" => {
                            circuit.add_behavioral_voltage_source(t[0].clone(), p, n, evaluator)?
                        }
                        "I" => {
                            circuit.add_behavioral_current_source(t[0].clone(), p, n, evaluator)?
                        }
                        _ => return Err(invalid("behavioral source output must be V or I")),
                    }
                }
                b'E' | b'G' if t.len() >= 4 => {
                    let (p, n) = (circuit.node(&t[1]), circuit.node(&t[2]));
                    let mode = upper(&t[3]);
                    let evaluator: Option<BehavioralEvaluator> = if mode.starts_with("POLY(") {
                        let close = mode.find(')').ok_or_else(|| {
                            invalid("malformed dependent-source POLY() declaration")
                        })?;
                        let digits: String = mode[5..close]
                            .trim_start()
                            .chars()
                            .take_while(char::is_ascii_digit)
                            .collect();
                        let dims: usize = digits
                            .parse()
                            .map_err(|_| invalid("stoul: no conversion"))?;
                        require(
                            dims == 1,
                            "dependent-source POLY baseline currently supports POLY(1)",
                        )?;
                        require(
                            t.len() >= 7,
                            "POLY(1) source requires control pair and coefficients",
                        )?;
                        let (cp, cn) = (circuit.node(&t[4]), circuit.node(&t[5]));
                        let coefficients = t[6..]
                            .iter()
                            .map(|s| value(s))
                            .collect::<Result<Vec<f64>>>()?;
                        Some(Arc::new(move |nodes: &[f64], _t: f64| {
                            let (Some(a), Some(b)) = (nodes.get(cp), nodes.get(cn)) else {
                                return Err(invalid("POLY control node out of range"));
                            };
                            let x = a - b;
                            Ok(coefficients.iter().rev().fold(0.0, |r, c| r * x + c))
                        }))
                    } else if mode == "TABLE" {
                        let spec = joined(&t, 4);
                        let eq = spec.find('=').ok_or_else(|| {
                            invalid("TABLE source requires input-expression = (x,y) points")
                        })?;
                        let input =
                            behavioral_evaluator(&mut circuit, trim(&spec[..eq]), &parameters)?;
                        let point_text: String = spec[eq + 1..]
                            .chars()
                            .map(|c| if matches!(c, '(' | ')' | ',') { ' ' } else { c })
                            .collect();
                        let pt = tokens(&point_text);
                        require(
                            pt.len() >= 4 && pt.len() % 2 == 0,
                            "TABLE source requires at least two (x,y) points",
                        )?;
                        let mut points = Vec::with_capacity(pt.len() / 2);
                        for pair in pt.chunks(2) {
                            points.push((value(&pair[0])?, value(&pair[1])?));
                        }
                        points.sort_by(|a, b| a.0.total_cmp(&b.0));
                        require(
                            points.windows(2).all(|w| w[1].0 > w[0].0),
                            "TABLE abscissas must be strictly increasing",
                        )?;
                        Some(Arc::new(move |nodes: &[f64], time: f64| {
                            let x = input(nodes, time)?;
                            let (first, last) = (points[0], points[points.len() - 1]);
                            if x <= first.0 {
                                return Ok(first.1);
                            }
                            if x >= last.0 {
                                return Ok(last.1);
                            }
                            let up = points.partition_point(|p| p.0 <= x);
                            let (lo, hi) = (points[up - 1], points[up]);
                            let f = (x - lo.0) / (hi.0 - lo.0);
                            Ok(lo.1 + f * (hi.1 - lo.1))
                        }))
                    } else {
                        None
                    };
                    match evaluator {
                        Some(e) if kind == b'E' => {
                            circuit.add_behavioral_voltage_source(t[0].clone(), p, n, e)?
                        }
                        Some(e) => circuit.add_behavioral_current_source(t[0].clone(), p, n, e)?,
                        None => {
                            require(
                                t.len() >= 6,
                                "dependent source requires control nodes and gain",
                            )?;
                            let (cp, cn) = (circuit.node(&t[3]), circuit.node(&t[4]));
                            let gain = value(&t[5])?;
                            if kind == b'G' {
                                circuit.add_vccs(t[0].clone(), p, n, cp, cn, gain)?;
                            } else {
                                circuit.add_vcvs(t[0].clone(), p, n, cp, cn, gain)?;
                            }
                        }
                    }
                }
                b'F' | b'H' if t.len() >= 5 => {
                    let (p, n) = (circuit.node(&t[1]), circuit.node(&t[2]));
                    let gain = value(&t[4])?;
                    if kind == b'F' {
                        circuit.add_cccs(t[0].clone(), p, n, &t[3], gain)?;
                    } else {
                        circuit.add_ccvs(t[0].clone(), p, n, &t[3], gain)?;
                    }
                }
                b'D' if t.len() >= 4 => {
                    let (a, k) = (circuit.node(&t[1]), circuit.node(&t[2]));
                    circuit.add_diode(t[0].clone(), a, k, &t[3])?;
                }
                b'M' if t.len() >= 6 => {
                    let n: Vec<usize> = (1..5).map(|i| circuit.node(&t[i])).collect();
                    circuit.add_mosfet(t[0].clone(), n[0], n[1], n[2], n[3], &t[5])?;
                }
                b'Q' if t.len() >= 5 => {
                    let n: Vec<usize> = (1..4).map(|i| circuit.node(&t[i])).collect();
                    circuit.add_bjt(t[0].clone(), n[0], n[1], n[2], &t[t.len() - 1])?;
                }
                b'J' if t.len() >= 5 => {
                    let n: Vec<usize> = (1..4).map(|i| circuit.node(&t[i])).collect();
                    circuit.add_jfet(t[0].clone(), n[0], n[1], n[2], &t[4])?;
                }
                b'S' if t.len() >= 6 => {
                    let n: Vec<usize> = (1..5).map(|i| circuit.node(&t[i])).collect();
                    circuit.add_switch(t[0].clone(), n[0], n[1], n[2], n[3], &t[5])?;
                }
                _ => {}
            }
        }
        for raw in &mutuals {
            let t = tokens(raw);
            if t.len() >= 4 {
                circuit.add_mutual_inductance(t[0].clone(), &t[1], &t[2], value(&t[3])?)?;
            }
        }
        for raw in &ics {
            apply_voltage_directive(&mut circuit, raw, true, &value)?;
        }
        for raw in &nodesets {
            apply_voltage_directive(&mut circuit, raw, false, &value)?;
        }
        Ok(circuit)
    }

    /// Parse a netlist file with recursive `.INCLUDE` and `.LIB <file> <section>` support;
    /// a non-empty `library_section` selects that `.LIB` section of the top file.
    pub fn parse_spice_file(path: impl AsRef<Path>, library_section: &str) -> Result<Circuit> {
        let mut active = HashSet::new();
        let lines =
            load_spice_file_recursive(&PathBuf::from(path.as_ref()), library_section, &mut active)?;
        let mut merged = String::new();
        for l in lines {
            merged.push_str(&l);
            merged.push('\n');
        }
        Circuit::parse_spice(&merged)
    }
}
