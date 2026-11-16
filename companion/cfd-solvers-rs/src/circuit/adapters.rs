//! External-model seams: the OSDI descriptor stamping adapter, OpenVAF and ngspice/Xyce command
//! plans, the SPICE compatibility profile, and OSDI module discovery. Ported from upstream
//! `cfd/circuit/{osdi_adapter, osdi_loader, external_spice_adapter, spice_compat}` (MIT,
//! cfd_solvers contributors).
//!
//! Upstream's `OsdiLibrary` `dlopen`s the module. Loading foreign code needs `unsafe`, which
//! this crate forbids, so [`OsdiLibrary::open`] instead reads the exported
//! `OSDI_VERSION_MAJOR/MINOR`, `OSDI_NUM_DESCRIPTORS` and `OSDI_DESCRIPTORS` symbols straight
//! from an ELF64 little-endian shared object. Descriptor evaluation stays behind
//! [`OsdiDescriptorAdapter`], as upstream.
use crate::{require, Result, SolverError};
use std::{path::Path, sync::Arc};

/// Fills terminal currents and a row-major Jacobian from terminal voltages.
pub type OsdiStampEvaluator =
    Arc<dyn Fn(&[f64], &mut [f64], &mut [f64]) -> Result<()> + Send + Sync>;

#[derive(Clone)]
pub struct OsdiDescriptorAdapter {
    terminal_count: usize,
    evaluator: OsdiStampEvaluator,
}

impl std::fmt::Debug for OsdiDescriptorAdapter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OsdiDescriptorAdapter")
            .field("terminal_count", &self.terminal_count)
            .finish_non_exhaustive()
    }
}

impl OsdiDescriptorAdapter {
    pub fn new(terminal_count: usize, evaluator: OsdiStampEvaluator) -> Result<Self> {
        require(terminal_count > 0, "invalid OSDI descriptor adapter")?;
        Ok(Self {
            terminal_count,
            evaluator,
        })
    }
    pub fn stamp(
        &self,
        terminal_voltage: &[f64],
        terminal_current: &mut [f64],
        jacobian: &mut [f64],
    ) -> Result<()> {
        let n = self.terminal_count;
        require(
            terminal_voltage.len() == n && terminal_current.len() == n && jacobian.len() == n * n,
            "OSDI descriptor stamp buffer size mismatch",
        )?;
        (self.evaluator)(terminal_voltage, terminal_current, jacobian)
    }
    pub fn terminal_count(&self) -> usize {
        self.terminal_count
    }
    /// The adapter as a circuit static-device evaluator.
    pub fn static_evaluator(&self) -> crate::tcad::StaticDeviceEvaluator {
        let adapter = self.clone();
        Arc::new(move |v: &[f64]| {
            let n = adapter.terminal_count;
            let mut current = vec![0.0; n];
            let mut jacobian = vec![0.0; n * n];
            adapter.stamp(v, &mut current, &mut jacobian)?;
            Ok(crate::tcad::StaticDeviceEvaluation {
                terminal_current: current,
                jacobian,
            })
        })
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct OpenVafCompilePlan {
    pub argv: Vec<String>,
    pub source: String,
    pub output_library: String,
}

/// `compiler --output <library> <source>` (upstream default compiler: `openvaf`).
pub fn make_openvaf_compile_plan(
    source: &str,
    output_library: &str,
    compiler: &str,
) -> Result<OpenVafCompilePlan> {
    require(
        !source.is_empty() && !output_library.is_empty() && !compiler.is_empty(),
        "OpenVAF compile plan paths must not be empty",
    )?;
    Ok(OpenVafCompilePlan {
        argv: vec![
            compiler.into(),
            "--output".into(),
            output_library.into(),
            source.into(),
        ],
        source: source.into(),
        output_library: output_library.into(),
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ExternalSpiceEngine {
    #[default]
    Ngspice,
    Xyce,
}

/// Command description for optional external cross-validation; never executed here.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ExternalSpiceRunPlan {
    pub engine: ExternalSpiceEngine,
    pub argv: Vec<String>,
    pub netlist: String,
    pub raw_output: String,
}

/// ngspice: `ngspice -b -r <raw> <netlist>` (extra arguments ignored, as upstream);
/// Xyce: `Xyce <extra...> <netlist>`.
pub fn make_external_spice_run_plan(
    engine: ExternalSpiceEngine,
    netlist: &str,
    raw_output: &str,
    extra_arguments: &[String],
) -> Result<ExternalSpiceRunPlan> {
    require(
        !netlist.is_empty() && !raw_output.is_empty(),
        "external SPICE plan paths must not be empty",
    )?;
    let argv = match engine {
        ExternalSpiceEngine::Ngspice => vec![
            "ngspice".into(),
            "-b".into(),
            "-r".into(),
            raw_output.into(),
            netlist.into(),
        ],
        ExternalSpiceEngine::Xyce => {
            let mut a = vec!["Xyce".to_string()];
            a.extend(extra_arguments.iter().cloned());
            a.push(netlist.into());
            a
        }
    };
    Ok(ExternalSpiceRunPlan {
        engine,
        argv,
        netlist: netlist.into(),
        raw_output: raw_output.into(),
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum SpiceCompatibilityMode {
    #[default]
    Native,
    Pspice,
    Ltspice,
    Hspice,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SpiceCompatibilityProfile {
    pub mode: SpiceCompatibilityMode,
    pub engineering_suffixes: bool,
    pub common_model_cards: bool,
    pub subcircuits_and_parameters: bool,
    pub supported_directives: Vec<&'static str>,
}

/// The parser's portable interchange subset (identical for every mode, as upstream).
pub fn spice_compatibility_profile(mode: SpiceCompatibilityMode) -> SpiceCompatibilityProfile {
    SpiceCompatibilityProfile {
        mode,
        engineering_suffixes: true,
        common_model_cards: true,
        subcircuits_and_parameters: true,
        supported_directives: vec![
            ".MODEL", ".SUBCKT", ".PARAM", ".FUNC", ".INCLUDE", ".LIB", ".IC", ".NODESET", ".TRAN",
            ".AC", ".DC", ".OP",
        ],
    }
}

/// ABI discovery for an OSDI module (OpenVAF output).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OsdiLibrary {
    path: String,
    version_major: u32,
    version_minor: u32,
    descriptor_count: u32,
    descriptors_file_offset: Option<u64>,
}

fn rd<const N: usize>(b: &[u8], at: u64) -> Option<[u8; N]> {
    let start = usize::try_from(at).ok()?;
    b.get(start..start.checked_add(N)?)?.try_into().ok()
}
fn u16_at(b: &[u8], at: u64) -> Option<u16> {
    rd::<2>(b, at).map(u16::from_le_bytes)
}
fn u32_at(b: &[u8], at: u64) -> Option<u32> {
    rd::<4>(b, at).map(u32::from_le_bytes)
}
fn u64_at(b: &[u8], at: u64) -> Option<u64> {
    rd::<8>(b, at).map(u64::from_le_bytes)
}

#[derive(Clone, Copy)]
struct Section {
    kind: u32,
    addr: u64,
    offset: u64,
    size: u64,
    link: u32,
    entsize: u64,
}

/// Exported symbol → (file offset of its storage or `None` for NOBITS, size).
fn elf_symbols(b: &[u8]) -> Option<Vec<(String, Option<u64>, u64)>> {
    if b.get(..4)? != b"\x7fELF" || *b.get(4)? != 2 || *b.get(5)? != 1 {
        return None;
    }
    let shoff = u64_at(b, 0x28)?;
    let shentsize = u64::from(u16_at(b, 0x3a)?);
    let shnum = u64::from(u16_at(b, 0x3c)?);
    if shentsize < 64 {
        return None;
    }
    let mut sections = Vec::new();
    for i in 0..shnum {
        let h = shoff.checked_add(i.checked_mul(shentsize)?)?;
        sections.push(Section {
            kind: u32_at(b, h + 4)?,
            addr: u64_at(b, h + 16)?,
            offset: u64_at(b, h + 24)?,
            size: u64_at(b, h + 32)?,
            link: u32_at(b, h + 40)?,
            entsize: u64_at(b, h + 56)?,
        });
    }
    let table = sections
        .iter()
        .find(|s| s.kind == 11)
        .or_else(|| sections.iter().find(|s| s.kind == 2))?;
    let strtab = sections.get(usize::try_from(table.link).ok()?)?;
    let entsize = if table.entsize == 0 {
        24
    } else {
        table.entsize
    };
    let mut out = Vec::new();
    for k in 0..table.size / entsize {
        let e = table.offset.checked_add(k * entsize)?;
        let name_off = u64::from(u32_at(b, e)?);
        let shndx = u16_at(b, e + 6)?;
        let value = u64_at(b, e + 8)?;
        let size = u64_at(b, e + 16)?;
        let start = usize::try_from(strtab.offset.checked_add(name_off)?).ok()?;
        let name_bytes = b.get(start..)?;
        let len = name_bytes.iter().position(|&c| c == 0)?;
        let name = String::from_utf8_lossy(&name_bytes[..len]).into_owned();
        if shndx == 0 || shndx >= 0xff00 {
            continue;
        }
        let sec = sections.get(usize::from(shndx))?;
        let file = if sec.kind == 8 {
            None
        } else {
            value
                .checked_sub(sec.addr)
                .and_then(|d| sec.offset.checked_add(d))
        };
        out.push((name, file, size));
    }
    Some(out)
}

impl OsdiLibrary {
    /// Read the OSDI ABI symbols of an ELF64 little-endian module without loading it.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path_ref = path.as_ref();
        let path_text = path_ref.to_string_lossy().into_owned();
        require(!path_text.is_empty(), "OSDI module path must not be empty")?;
        let bytes = std::fs::read(path_ref)
            .map_err(|e| SolverError::Io(format!("cannot load OSDI module {path_text}: {e}")))?;
        Self::from_bytes(&path_text, &bytes)
    }

    /// Same as [`OsdiLibrary::open`] for an in-memory image.
    pub fn from_bytes(path: &str, bytes: &[u8]) -> Result<Self> {
        let symbols = elf_symbols(bytes).ok_or_else(|| {
            SolverError::Unsupported(format!(
                "cannot load OSDI module {path}: not a readable ELF64 little-endian shared object"
            ))
        })?;
        let find = |name: &str| symbols.iter().find(|s| s.0 == name);
        let required_u32 = |name: &str| -> Result<u32> {
            let (_, offset, _) = find(name).ok_or_else(|| {
                SolverError::InvalidInput(format!("OSDI module missing symbol {name}"))
            })?;
            Ok(match offset {
                None => 0,
                Some(o) => u32_at(bytes, *o).ok_or_else(|| {
                    SolverError::InvalidInput(format!("OSDI symbol {name} is outside the file"))
                })?,
            })
        };
        let version_major = required_u32("OSDI_VERSION_MAJOR")?;
        let version_minor = required_u32("OSDI_VERSION_MINOR")?;
        let descriptor_count = required_u32("OSDI_NUM_DESCRIPTORS")?;
        let descriptors = find("OSDI_DESCRIPTORS").map(|s| s.1.unwrap_or(0));
        require(
            descriptor_count == 0 || descriptors.is_some(),
            "OSDI module declares descriptors but exports no OSDI_DESCRIPTORS symbol",
        )?;
        Ok(Self {
            path: path.into(),
            version_major,
            version_minor,
            descriptor_count,
            descriptors_file_offset: descriptors,
        })
    }
    pub fn version_major(&self) -> u32 {
        self.version_major
    }
    pub fn version_minor(&self) -> u32 {
        self.version_minor
    }
    pub fn descriptor_count(&self) -> u32 {
        self.descriptor_count
    }
    /// File offset of the `OSDI_DESCRIPTORS` table (upstream exposes the loaded pointer).
    pub fn descriptors_file_offset(&self) -> Option<u64> {
        self.descriptors_file_offset
    }
    pub fn path(&self) -> &str {
        &self.path
    }
}
