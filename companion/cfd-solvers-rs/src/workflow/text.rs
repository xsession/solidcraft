//! Small text helpers reproducing the C/C++ library behaviour the upstream workflow code relies on
//! (`isspace`, `std::getline`, `std::stod`, `std::stoull`, `std::stoi`, iostream `%.Ng`).
use crate::{Result, SolverError};

/// C-locale `isspace`: space, `\t`, `\n`, `\v`, `\f`, `\r`.
pub(crate) fn c_isspace(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

pub(crate) fn format_g17(v: f64) -> String {
    crate::lbm::format_general(v, 17)
}

/// Pieces as produced by repeated `std::getline(stream, piece, sep)`: a trailing empty piece
/// (text ending with the separator, or empty text) is not produced.
pub(crate) fn getline_split(text: &str, sep: u8) -> Vec<&str> {
    let mut out: Vec<&str> = text.split(sep as char).collect();
    if out.last().is_some_and(|s| s.is_empty()) {
        out.pop();
    }
    out
}

/// Whitespace-separated tokens (`stream >> token`).
pub(crate) fn ws_tokens(text: &str) -> Vec<String> {
    text.split(|c: char| c.is_ascii() && c_isspace(c as u8))
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect()
}

pub(crate) fn trim_ascii(s: &str) -> &str {
    s.trim_matches(|c| matches!(c, ' ' | '\t' | '\r' | '\n'))
}

/// `std::stod` on a whole token: leading whitespace allowed, the whole remainder must parse.
/// Hexadecimal floats are not accepted.
pub(crate) fn stod_full(token: &str) -> Option<f64> {
    let t = token.trim_start_matches(|c: char| c.is_ascii() && c_isspace(c as u8));
    if t.is_empty() || t.starts_with("0x") || t.starts_with("0X") {
        return None;
    }
    t.parse::<f64>().ok()
}

/// Longest numeric prefix as read by `std::stod`/`strtod`; no prefix is an error.
pub(crate) fn strtod_prefix(token: &str) -> Result<f64> {
    let t = token.trim_start_matches(|c: char| c.is_ascii() && c_isspace(c as u8));
    for end in (1..=t.len()).rev() {
        if let Some(prefix) = t.get(..end) {
            if let Ok(v) = prefix.parse::<f64>() {
                return Ok(v);
            }
        }
    }
    Err(SolverError::InvalidInput(format!(
        "stod: no conversion for '{token}'"
    )))
}

/// `std::stoull` on a whole token (leading whitespace and sign allowed; a minus sign wraps
/// modulo 2^64 as in C). Overflow is a failure.
pub(crate) fn stoull_full(token: &str) -> Option<u64> {
    let t = token.trim_start_matches(|c: char| c.is_ascii() && c_isspace(c as u8));
    let (neg, digits) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
        return None;
    }
    let v: u64 = digits.parse().ok()?;
    Some(if neg { 0u64.wrapping_sub(v) } else { v })
}

pub(crate) fn parse_size(token: &str) -> Option<usize> {
    stoull_full(token).and_then(|v| usize::try_from(v).ok())
}

/// `std::stoi` (numeric prefix, leading whitespace and sign); `None` when it would throw.
pub(crate) fn stoi_prefix(text: &str) -> Option<i32> {
    let t = text.trim_start_matches(|c: char| c.is_ascii() && c_isspace(c as u8));
    let b = t.as_bytes();
    let mut end = usize::from(matches!(b.first(), Some(b'-' | b'+')));
    let digits_start = end;
    while end < b.len() && b[end].is_ascii_digit() {
        end += 1;
    }
    if end == digits_start {
        return None;
    }
    t[..end].parse::<i32>().ok()
}

pub(crate) fn upper(s: &str) -> String {
    s.to_ascii_uppercase()
}

pub(crate) fn lower(s: &str) -> String {
    s.to_ascii_lowercase()
}

/// Message text of an error without the `kind:` prefix, i.e. what upstream's `e.what()` returns.
pub(crate) fn error_text(e: &SolverError) -> String {
    match e {
        SolverError::InvalidInput(s)
        | SolverError::Singular(s)
        | SolverError::NumericalFailure(s)
        | SolverError::Unsupported(s)
        | SolverError::Io(s) => s.clone(),
    }
}

pub(crate) fn io_err(what: impl std::fmt::Display) -> impl FnOnce(std::io::Error) -> SolverError {
    move |e| SolverError::Io(format!("{what}: {e}"))
}

/// Path text with `/` separators (`generic_string()`).
pub(crate) fn path_text(p: &std::path::Path) -> String {
    let s = p.to_string_lossy();
    if std::path::MAIN_SEPARATOR == '/' {
        s.into_owned()
    } else {
        s.replace(std::path::MAIN_SEPARATOR, "/")
    }
}

/// Add the owner-execute bit (no-op on non-Unix platforms).
pub(crate) fn add_owner_exec(path: &std::path::Path) -> Result<()> {
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let meta = std::fs::metadata(path).map_err(io_err("cannot stat file"))?;
        let mut perm = meta.permissions();
        perm.set_mode(perm.mode() | 0o100);
        std::fs::set_permissions(path, perm).map_err(io_err("cannot set permissions"))?;
    }
    #[cfg(not(unix))]
    let _ = path;
    Ok(())
}

pub(crate) fn write_file(path: &std::path::Path, text: &str, what: &str) -> Result<()> {
    std::fs::write(path, text).map_err(|e| SolverError::Io(format!("{what}: {e}")))
}

pub(crate) fn create_dirs(path: &std::path::Path) -> Result<()> {
    if path.as_os_str().is_empty() {
        return Ok(());
    }
    std::fs::create_dir_all(path).map_err(io_err(format!("cannot create {}", path.display())))
}
