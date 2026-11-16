// Test oracle only. Exercises upstream workflow/{schema, provenance, validation_catalog, campaign,
// campaign_control, deploy}. File-writing functions write under argv[1].
#include "cfd/workflow/campaign.hpp"
#include "cfd/workflow/campaign_control.hpp"
#include "cfd/workflow/deploy.hpp"
#include "cfd/workflow/provenance.hpp"
#include "cfd/workflow/schema.hpp"
#include "cfd/workflow/validation_catalog.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace cfd::workflow;
namespace fs = std::filesystem;
static bool first = true;
static std::string esc(const std::string& s) {
    std::ostringstream o;
    for (unsigned char c : s) {
        if (c == '"' || c == '\\') o << '\\' << c;
        else if (c == '\n') o << "\\n";
        else if (c == '\t') o << "\\t";
        else if (c == '\r') o << "\\r";
        else if (c < 0x20) o << "\\u" << std::hex << std::setw(4) << std::setfill('0') << (int)c << std::dec;
        else o << c;
    }
    return o.str();
}
static void key(const std::string& k) { std::cout << (first ? "" : ",\n") << "\"" << k << "\":"; first = false; }
static void str(const std::string& k, const std::string& v) { key(k); std::cout << "\"" << esc(v) << "\""; }
static void strs(const std::string& k, const std::vector<std::string>& v) {
    key(k); std::cout << "[";
    for (std::size_t i = 0; i < v.size(); ++i) std::cout << (i ? "," : "") << "\"" << esc(v[i]) << "\"";
    std::cout << "]";
}
static void vec(const std::string& k, const std::vector<double>& v) {
    key(k); std::cout << "[";
    for (std::size_t i = 0; i < v.size(); ++i) {
        std::cout << (i ? "," : "");
        if (std::isfinite(v[i])) std::cout << std::setprecision(17) << v[i];
        else std::cout << (std::isnan(v[i]) ? "\"nan\"" : (v[i] > 0 ? "\"inf\"" : "\"-inf\""));
    }
    std::cout << "]";
}
static std::string slurp(const fs::path& p) { std::ifstream in(p, std::ios::binary); std::ostringstream o; o << in.rdbuf(); return o.str(); }
static std::string err(auto f) { try { f(); return ""; } catch (const std::exception& e) { return e.what(); } }
static std::string summary_text(const CampaignSummary& s) {
    std::ostringstream o; o << std::setprecision(17) << s.pending << ' ' << s.running << ' ' << s.done << ' ' << s.failed << ' ' << s.stopped << ' ' << s.best_objective << ' ' << s.best_case_id; return o.str();
}
static std::vector<std::string> cases_text(const std::vector<CampaignCase>& cs) {
    std::vector<std::string> out;
    for (auto& c : cs) {
        std::vector<std::string> keys; for (auto& kv : c.values) keys.push_back(kv.first); std::sort(keys.begin(), keys.end());
        std::string s = c.case_id; for (auto& k : keys) s += "|" + k + "=" + c.values.at(k); out.push_back(s);
    }
    return out;
}
static std::vector<std::string> registry_text(const std::vector<CampaignRegistryEntry>& r) {
    std::vector<std::string> out;
    for (auto& e : r) { std::ostringstream o; o << std::setprecision(17) << e.case_id << '|' << campaign_status_name(e.status) << '|' << e.objective << '|' << e.iterations << '|' << e.message; out.push_back(o.str()); }
    return out;
}
static std::vector<std::string> cmds_text(const std::vector<RemoteCommandPlan>& cs) {
    std::vector<std::string> out;
    for (auto& c : cs) { std::string s = c.server_name + "|" + c.description + "|" + (c.uses_ssh ? "ssh" : "local") + "|" + c.display_command; for (auto& a : c.argv) s += "|" + a; out.push_back(s); }
    return out;
}
int main(int argc, char** argv) {
    const fs::path root = argc > 1 ? fs::path(argv[1]) : fs::path("/tmp/cfdwf_oracle");
    fs::remove_all(root);
    std::cout << "{\n";
    // schema
    CaseConfig cfg;
    cfg.set("velocity", {1.25, "m/s"}); cfg.set("density", {998.2, "kg/m^3"}); cfg.set("a\"b", {1e-7, "x\\y"}); cfg.set("temperature", {300.0 / 7.0, "K"});
    str("case_json", cfg.to_json());
    str("validate_ok", err([&] { cfg.validate({{"velocity", "m/s", 0.0, 10.0, true}, {"missing", "", {}, {}, false}}); }));
    str("validate_missing", err([&] { cfg.validate({{"missing", "", {}, {}, true}}); }));
    str("validate_unit", err([&] { cfg.validate({{"velocity", "km/h", {}, {}, true}}); }));
    str("validate_min", err([&] { cfg.validate({{"velocity", "", 2.0, {}, true}}); }));
    str("validate_max", err([&] { cfg.validate({{"velocity", "", {}, 1.0, true}}); }));
    auto parsed = case_config_from_json("{\"parameters\":{\"p\":{\"value\":101325,\"unit\":\"Pa\"}, \"q\" : { \"value\" : -1.5e-3 , \"unit\" : \"\" }, \"bad\":{\"value\":\"x\",\"unit\":\"m\"}, \"r\":{\"value\":2.5e,\"unit\":\"s\"}}}");
    str("from_json", parsed.to_json());
    auto yaml = case_config_from_yaml("# comment\nvelocity: 2.5 m/s\n  pressure : -1e5 Pa  \nbad: 3\nunit_hash: 4 m#x\nT: 300 K\r\nlast:7 s");
    str("from_yaml", yaml.to_json());
    MaterialDatabase db;
    db.add({"water", {{"rho", {998.2, "kg/m^3"}}, {"k", {0.6, "W/m/K"}}}});
    db.add({"air", {{"rho", {1.225, "kg/m^3"}}}});
    db.add({"empty", {}});
    str("materials", db.to_json());
    str("material_empty_name", err([&] { db.add({"", {}}); }));
    ResultProvenance prov{"fvm", "0.4.0", "abc123", "2026-11-14T12:34:56Z", {{"zeta", "z"}, {"alpha", "q\"uote"}}};
    str("provenance", provenance_json(prov));
    { std::vector<std::string> v; for (auto& c : validation_catalog()) v.push_back(c.name + "|" + c.family + "|" + c.command + "|" + c.criterion); strs("catalog", v); }
    // campaign generation
    CampaignRandom rng; std::vector<double> u; for (int i = 0; i < 8; ++i) u.push_back(rng.uniform_open01()); vec("uniform", u);
    std::vector<CampaignParameter> disc{{"mesh", CampaignParameterKind::discrete, 0, 0, {"coarse", "fine"}}, {"model", CampaignParameterKind::discrete, 0, 0, {"kEps", "kOmega", "SA"}}};
    strs("factorial", cases_text(generate_factorial_campaign(disc)));
    std::vector<CampaignParameter> mixed{{"u", CampaignParameterKind::continuous, 0.5, 2.5, {}}, {"T", CampaignParameterKind::continuous, 280, 320, {}}, {"mesh", CampaignParameterKind::discrete, 0, 0, {"a", "b", "c"}}};
    strs("lhs", cases_text(generate_latin_hypercube_campaign(mixed, 7)));
    strs("lhs_seed", cases_text(generate_latin_hypercube_campaign(mixed, 5, CampaignRandom{12345})));
    str("factorial_err", err([&] { (void)generate_factorial_campaign(mixed); }));
    str("reserved_err", err([&] { std::vector<CampaignParameter> p{{"case_id", CampaignParameterKind::discrete, 0, 0, {"x"}}}; (void)generate_factorial_campaign(p); }));
    CampaignCase row{"case0001", {{"mesh", "fine"}, {"u", "1.5"}, {"turb", "1"}, {"off", "false"}}};
    std::string templ = "u={u} mesh={mesh} {unknown_brace\n<!-- IF turb -->\nturbulence on\n<!-- IF mesh=fine -->\nfine {mesh}{u}\n<!-- ENDIF -->\n<!-- ENDIF -->\n<!-- IF off -->\nhidden\n<!-- ENDIF -->\n<!-- IF mesh != coarse -->\nnot coarse\n<!-- ENDIF -->\n{{u}}\nend";
    str("render", render_campaign_template(templ, row, false));
    str("render_strict_err", err([&] { (void)render_campaign_template("x {nope}", row, true); }));
    str("render_unclosed", err([&] { (void)render_campaign_template("<!-- IF u -->\nx", row, true); }));
    str("render_endif", err([&] { (void)render_campaign_template("<!-- ENDIF -->", row, true); }));
    // registry
    std::vector<CampaignRegistryEntry> reg;
    update_case_status(reg, {"case0001", CampaignCaseStatus::done, 0.25, 12, "ok"});
    update_case_status(reg, {"case0002", CampaignCaseStatus::done, 0.125, 7, "tab\there\nnew\\line"});
    update_case_status(reg, {"case0003", CampaignCaseStatus::failed, 1.0 / 3.0, 0, ""});
    update_case_status(reg, {"case0001", CampaignCaseStatus::done, 0.5, 13, "again"});
    update_case_status(reg, {"case0004", CampaignCaseStatus::running, 0, 0, "r"});
    str("summary", summary_text(summarize_campaign_registry(reg)));
    save_campaign_registry(root / "reg" / "registry.tsv", reg);
    str("registry_tsv", slurp(root / "reg" / "registry.tsv"));
    strs("registry_loaded", registry_text(load_campaign_registry(root / "reg" / "registry.tsv")));
    { std::vector<std::string> v; for (auto n : {"pending", "RUNNING", "Success", "completed", "error", "cancelled", "STOPPED"}) v.push_back(campaign_status_name(campaign_status_from_name(n))); strs("status_names", v); }
    // case folders
    std::vector<CampaignTemplateFile> tf{{"system/controlDict", "u {u}\nmesh {mesh}\n"}, {"run.sh", "#!/bin/sh\necho {mesh}\n<!-- IF mesh=b -->\necho b\n<!-- ENDIF -->\n"}};
    auto lcases = generate_latin_hypercube_campaign(mixed, 3);
    auto pr = write_campaign_case_folders(root / "camp", lcases, tf);
    vec("persist", {(double)pr.cases_written, (double)pr.files_written, (double)pr.registry_entries});
    str("persist_registry", slurp(root / "camp" / "registry.tsv"));
    str("persist_manifest", slurp(root / "camp" / "template_manifest.tsv"));
    str("persist_doe", slurp(root / "camp" / "case0002" / "doe_row.csv"));
    str("persist_ctrl", slurp(root / "camp" / "case0002" / "system" / "controlDict"));
    str("persist_run", slurp(root / "camp" / "case0003" / "run.sh"));
    str("persist_exists_err", err([&] { (void)write_campaign_case_folders(root / "camp", lcases, tf); }).substr(0, 37));
    str("persist_escape_err", err([&] { std::vector<CampaignTemplateFile> b{{"../x", ""}}; (void)write_campaign_case_folders(root / "camp2", lcases, b); }));
    // command plans
    SolverAdapterDescriptor ad{"cfd-solve", "cfd-solve", "/opt/cfd/bin/cfd-solve", {SolverCapability::residuals, SolverCapability::control}, {"stop", "checkpoint"}};
    std::vector<std::string> plans;
    for (auto rt : {SolverRuntime::native, SolverRuntime::docker, SolverRuntime::singularity, SolverRuntime::slurm}) {
        SolverRunRequest req{ad, "/data/my case", rt, 4, 2, {"--flag", "it's"}, "img:1", "img.sif", "gpu"};
        auto p = build_solver_command_plan(req);
        std::string s = std::to_string((int)p.runtime) + "|" + (p.requires_scheduler ? "1" : "0") + "|" + p.display_command;
        for (auto& [k, v] : p.environment) s += "|" + k + "=" + v;
        for (auto& a : p.argv) s += "|" + a;
        plans.push_back(s);
    }
    strs("plans", plans);
    str("plan_docker_err", err([&] { SolverRunRequest req{ad, "/x", SolverRuntime::docker}; (void)build_solver_command_plan(req); }));
    // log parsing
    std::string log = "iter 1 res 1e-1\nIteration=2, residual=5e-2 field=U\n3 0.025 extra\nEQUATION p ITER 4 RES 1.25e-2\nnoise 5\n7 inf\nperf.wall_time 12.5 s\nWALL_TIME: 3\nCELLS_PER_S=1e6 cells/s\nperf.mem 1.5\nother 1\n";
    { std::vector<double> v; for (auto& s : parse_residual_history(log)) { v.push_back(s.iteration); v.push_back(s.residual); } vec("residuals", v); std::vector<std::string> e; for (auto& s : parse_residual_history(log)) e.push_back(s.equation); strs("residual_eq", e); }
    { std::vector<std::string> v; for (auto& m : parse_performance_metrics(log)) { std::ostringstream o; o << std::setprecision(17) << m.key << '|' << m.value << '|' << m.units; v.push_back(o.str()); } strs("perf", v); }
    { std::vector<std::string> v; for (auto l : {"all done", "Fatal: x", "job cancelled", "iter 3", "SUCCESS but error"}) v.push_back(campaign_status_name(detect_solver_outcome_from_log(l))); strs("outcomes", v); }
    // outputs discovery + refresh
    fs::create_directories(root / "camp" / "case0001" / "logs");
    std::ofstream(root / "camp" / "case0001" / "logs" / "solver.log") << "iter 1 res 0.5\niter 2 res 0.25\n";
    std::ofstream(root / "camp" / "case0001" / "residual_history.dat") << "3 0.125\n";
    std::ofstream(root / "camp" / "case0001" / "timing.txt") << "perf.wall 2 s\n";
    std::ofstream(root / "camp" / "case0001" / "notes.md") << "DONE\n";
    std::ofstream(root / "camp" / "case0002" / "run.out") << "iter 9 res 1e-3\nFAILED\n";
    { auto d = discover_campaign_outputs(root / "camp" / "case0001"); std::vector<std::string> v; for (auto& p : d.log_files) v.push_back("L" + fs::relative(p, root).generic_string()); for (auto& p : d.residual_files) v.push_back("R" + fs::relative(p, root).generic_string()); for (auto& p : d.performance_files) v.push_back("P" + fs::relative(p, root).generic_string()); strs("discover", v); }
    str("refresh_summary", summary_text(refresh_campaign_registry_from_outputs(root / "camp")));
    strs("refresh_registry", registry_text(load_campaign_registry(root / "camp" / "registry.tsv")));
    // control
    CampaignControlRequest creq{"case0002", CampaignControlAction::checkpoint, "100\nsteps", "user\\req"};
    auto dir = write_campaign_control_directive(root / "camp", ad, creq);
    str("directive", slurp(dir.path));
    str("directive_name", fs::relative(dir.path, root).generic_string());
    str("control_err", err([&] { (void)write_campaign_control_directive(root / "camp", ad, {"case0002", CampaignControlAction::extend, "", ""}); }));
    { auto ds = discover_campaign_control_directives(root / "camp" / "case0002"); std::vector<std::string> v; for (auto& d : ds) v.push_back(d.case_id + "|" + campaign_control_action_name(d.action) + "|" + d.value + "|" + d.reason + "|" + d.path.filename().string()); strs("directives", v); }
    auto jobs = parse_slurm_queue_table("JOBID NAME ST TIME\n# c\n101,case0001,R,12.5\n102 case0005 PD\n103\tcase0004\tCD\t7\n104 case0003 WEIRD x\nshort row\n");
    { std::vector<std::string> v; for (auto& j : jobs) { std::ostringstream o; o << std::setprecision(17) << j.job_id << '|' << j.case_id << '|' << scheduler_job_state_name(j.state) << '|' << j.elapsed_seconds << '|' << j.raw_state; v.push_back(o.str()); } strs("slurm", v); }
    { auto r2 = reg; str("sched_summary", summary_text(apply_scheduler_records_to_registry(r2, jobs))); strs("sched_registry", registry_text(r2)); }
    // optimisation
    auto rosen = [](std::span<const double> x) { return (1 - x[0]) * (1 - x[0]) + 10 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]); };
    vec("fd_central", finite_difference_gradient(rosen, std::vector<double>{0.3, -0.2}));
    vec("fd_forward", finite_difference_gradient(rosen, std::vector<double>{30.0, -0.2}, {1e-6, 1e-4, false}));
    { auto r = run_gradient_descent_campaign(rosen, {-0.5, 0.5}, {200, 0.05, 1e-8, {}}); auto v = r.parameters; v.push_back(r.objective); v.push_back(r.iterations); v.push_back(r.converged); v.insert(v.end(), r.last_gradient.begin(), r.last_gradient.end()); vec("gd", v); }
    { auto q = [](std::span<const double> x) { return (x[0] - 3) * (x[0] - 3) + 2 * (x[1] + 1) * (x[1] + 1); }; auto r = run_gradient_descent_campaign(q, {0, 0}); auto v = r.parameters; v.push_back(r.objective); v.push_back(r.iterations); v.push_back(r.converged); vec("gd_quad", v); }
    // deploy: multi-server
    std::vector<CampaignServerDescriptor> servers{
        {"gpu-a", "10.0.0.11", "cfd", "/srv/cfd/campaigns", CampaignServerRuntime::docker, 2, {"gpu", "linux"}, "cfd-solvers:latest", true},
        {"cpu-b", "localhost", "", "/scratch/my runs", CampaignServerRuntime::native, 1, {"linux"}, "", true},
        {"cpu-a", "cluster.example", "", "/srv/cfd", CampaignServerRuntime::native, 3, {"linux"}, "", true},
        {"off", "10.0.0.99", "", "/srv", CampaignServerRuntime::native, 4, {"linux"}, "", false}};
    std::vector<CampaignRegistryEntry> preg;
    for (int i = 1; i <= 8; ++i) { std::ostringstream id; id << "case" << std::setw(4) << std::setfill('0') << i; preg.push_back({id.str(), i == 3 ? CampaignCaseStatus::done : CampaignCaseStatus::pending, 0, 0, ""}); }
    MultiServerSchedulingOptions so; so.max_cases_per_server = 2;
    auto mplan = plan_multi_server_campaign("/campaigns/local", preg, ad, servers, so);
    { std::vector<std::string> v; for (auto& a : mplan.assignments) v.push_back(a.case_id + "|" + a.server_name + "|" + a.local_case_dir.string() + "|" + a.remote_case_dir.string() + "|" + a.solver_plan.display_command); strs("ms_assign", v); }
    strs("ms_cmds", cmds_text(mplan.commands));
    vec("ms_counts", {(double)mplan.unassigned_cases, (double)mplan.active_servers});
    { MultiServerSchedulingOptions t; t.required_tag = "gpu"; auto p = plan_multi_server_campaign("/c", preg, ad, servers, t); vec("ms_tag", {(double)p.assignments.size(), (double)p.unassigned_cases, (double)p.active_servers}); }
    str("ms_err", err([&] { MultiServerSchedulingOptions t; t.required_tag = "nope"; (void)plan_multi_server_campaign("/c", preg, ad, servers, t); }));
    auto pw = write_multiserver_plan_files(root / "ms", mplan);
    str("ms_assign_tsv", slurp(pw.assignments_path)); str("ms_cmd_sh", slurp(pw.commands_path));
    auto ex = plan_multiserver_execution(mplan, servers);
    { std::vector<std::string> v; for (auto& h : ex.health_checks) v.push_back(h.server_name + "|" + h.check_name + "|" + (h.required ? "1" : "0") + "|" + h.command.display_command); strs("ex_health", v); }
    { std::vector<RemoteCommandPlan> c; for (auto& j : ex.jobs) { c.push_back(j.launch_command); c.push_back(j.status_command); c.push_back(j.cancel_command); c.push_back(j.fetch_logs_command); } strs("ex_jobs", cmds_text(c)); }
    vec("ex_counts", {(double)ex.native_jobs, (double)ex.docker_jobs, (double)ex.remote_jobs});
    auto ew = write_multiserver_execution_files(root / "ms", ex);
    str("ex_health_sh", slurp(ew.health_script_path)); str("ex_launch_sh", slurp(ew.launch_script_path)); str("ex_status_sh", slurp(ew.status_script_path)); str("ex_fetch_sh", slurp(ew.fetch_logs_script_path)); str("ex_jobs_tsv", slurp(ew.jobs_path));
    auto rs = parse_remote_job_status_table("case_id\tserver\tstate\texit\ncase0001\tgpu-a\tdone\t0\ncase0002\tgpu-a\tfailed\t2x\ncase0004\tcpu-b\tmissing\ncase0009\tcpu-a\tqueued\tzz\n");
    { std::vector<std::string> v; for (auto& r : rs) v.push_back(r.case_id + "|" + r.server_name + "|" + remote_job_state_name(r.state) + "|" + std::to_string(r.exit_code) + "|" + r.raw_state); strs("remote_status", v); }
    { auto r2 = preg; str("remote_summary", summary_text(apply_remote_job_status_to_registry(r2, rs))); strs("remote_registry", registry_text(r2)); }
    str("remote_err", err([&] { (void)parse_remote_job_status_table("a\tb\tweird\n"); }));
    { std::vector<std::string> mk{"TOKEN=", "PASSWORD=", "SECRET=", "KEY=", "AWS_", "GITHUB_", ""}; strs("redact", {redact_sensitive_command_display("env TOKEN=abc123 PASSWORD='p w' AWS_SECRET=xyz;run GITHUB_TOKEN=t KEY=", mk), redact_sensitive_command_display("plain command", mk)}); }
    RemoteSupervisionOptions sup;
    auto ex2 = ex; ex2.jobs[0].launch_command.display_command += " TOKEN=supersecret";
    auto sp = plan_multiserver_supervision(ex2, servers, sup);
    { std::vector<std::string> v; for (auto& p : sp.access_probes) v.push_back(p.server_name + "|" + p.probe_name + "|" + p.command.display_command); strs("sup_probes", v); }
    { std::vector<RemoteCommandPlan> c; for (auto& r : sp.retry_launches) c.push_back(r.command); for (auto& t : sp.log_tails) { c.push_back(t.stdout_tail_command); c.push_back(t.stderr_tail_command); } strs("sup_cmds", cmds_text(c)); }
    vec("sup_counts", {(double)sp.sensitive_commands_redacted, (double)sp.dashboard_cases.size()});
    str("dashboard_json", build_multiserver_dashboard_json(sp.dashboard_cases));
    auto sw = write_multiserver_supervision_files(root / "ms", sp);
    str("sup_access_sh", slurp(sw.access_script_path)); str("sup_retry_sh", slurp(sw.retry_launch_script_path)); str("sup_tail_sh", slurp(sw.tail_logs_script_path)); str("sup_tsv", slurp(sw.summary_path));
    auto cases = sp.dashboard_cases; cases[1].state = "Running"; cases[2].state = "stopped"; cases[3].state = "weird\x01"; cases[0].exit_code = -3;
    str("ctl_status", build_multiserver_controller_status_json(cases));
    str("ctl_events", build_multiserver_controller_events_ndjson(cases));
    RemoteControllerConfig cc; cc.dashboard_title = "My \"Campaign\""; cc.control_actions = {"stop", "flush"}; cc.port = 9443;
    auto cp = plan_multiserver_controller(sp, cc);
    str("ctl_openapi", cp.openapi_json); str("ctl_script", cp.controller_script); str("ctl_html", cp.index_html); str("ctl_js", cp.dashboard_js); str("ctl_css", cp.dashboard_css);
    { RemoteControllerConfig c2; c2.enable_event_history = false; c2.allow_mutating_actions = false; c2.require_authentication = false; c2.enable_tls = false; auto p2 = plan_multiserver_controller(sp, c2); str("ctl_openapi2", p2.openapi_json); str("ctl_script2", p2.controller_script); }
    auto cw = write_multiserver_controller_files(root / "ms", cp);
    for (auto [k, p] : std::vector<std::pair<std::string, fs::path>>{{"routes", cw.routes_path}, {"env", cw.env_path}, {"readme", cw.readme_path}, {"tokens", cw.tokens_example_path}, {"systemd", cw.systemd_path}, {"caddy", cw.reverse_proxy_path}, {"start", cw.start_script_path}}) str("ctlw_" + k, slurp(p));
    str("ctl_err", err([&] { RemoteControllerConfig c3; c3.port = 0; (void)plan_multiserver_controller(sp, c3); }));
    // docker
    DockerDeployConfig dc; dc.worker_replicas = 3; dc.worker_cpus = 4; dc.worker_memory = "8g"; dc.enable_gpu_runtime = true; dc.environment = {{"OMP_NUM_THREADS", "4"}};
    { std::vector<std::string> v; for (auto& f : generate_docker_deploy_files(dc)) v.push_back(f.relative_path.generic_string() + "\x1f" + f.contents); strs("docker_full", v); }
    { DockerDeployConfig d2; d2.include_manager_service = false; std::vector<std::string> v; for (auto& f : generate_docker_deploy_files(d2)) v.push_back(f.relative_path.generic_string() + "\x1f" + f.contents); strs("docker_min", v); }
    { auto r = write_docker_deploy_system(root / "docker", dc); vec("docker_written", {(double)r.files_written}); }
    std::cout << "\n}\n";
}
