// Test oracle only. Never linked into the Rust library or executable.
// The implementations exercised here remain in the MIT-licensed upstream checkout.
#include "cfd/solvers/fem/poisson1d.hpp"
#include "cfd/solvers/lbm/d2q9.hpp"
#include "cfd/solvers/fdtd/maxwell1d.hpp"
#include "cfd/solvers/fvm/compressible1d.hpp"
#include "cfd/solvers/fvm/diffusion2d.hpp"
#include "cfd/solvers/fem/poisson2d.hpp"
#include "cfd/solvers/fem/poisson3d.hpp"
#include "cfd/solvers/fem/heat2d.hpp"
#include "cfd/solvers/fem/elasticity2d.hpp"
#include "cfd/solvers/fem/modal_bar1d.hpp"
#include "cfd/solvers/fem/piezoelectric_bar1d.hpp"
#include "cfd/fem/material_laws.hpp"
#include "cfd/rf/network.hpp"
#include "cfd/solvers/optics/sequential.hpp"
#include "cfd/solvers/optics/gaussian_beam.hpp"
#include "cfd/solvers/optics/materials.hpp"
#include "cfd/solvers/optics/polarization.hpp"
#include "cfd/chemistry/kinetics.hpp"
#include "cfd/chemistry/thermo.hpp"
#include "cfd/chemistry/equilibrium.hpp"
#include "cfd/chemistry/implicit_reactor.hpp"
#include "cfd/acoustics/kspace2d.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

bool first=true;
void emit(const std::string& key,const std::vector<double>& values){
    if(!first)std::cout<<',';first=false;std::cout<<'"'<<key<<"\":[";
    for(std::size_t i=0;i<values.size();++i){if(i)std::cout<<',';std::cout<<values[i];}std::cout<<']';
}
int main(){
    std::cout<<std::setprecision(17)<<'{';
    cfd::fem::Poisson1D poisson(19,1.7);poisson.solve([](double x){return 2.0+x*x;},0.2,0.7);
    emit("poisson1d",{poisson.u().begin(),poisson.u().end()});
    for(int scenario=0;scenario<4;++scenario){
        cfd::lbm::D2Q9Config cfg;cfg.nx=12;cfg.ny=10;cfg.tau=0.73F;
        if(scenario==1)cfg.acceleration_x=1e-5F;
        cfd::lbm::D2Q9Solver solver(cfg);
        if(scenario==0)solver.initialize_taylor_green(0.023F);
        else{for(std::size_t x=0;x<12;++x){solver.set_solid(x,0);solver.set_wall_velocity(x,9,0.017F,0.0F);}solver.initialize_uniform();}
        if(scenario==2){solver.set_velocity_inlet_left(0.01F);solver.set_pressure_outlet_right(1.0F);}
        if(scenario==3){for(std::size_t x=0;x<12;++x){solver.set_interpolated_wall_link(x,1,2,0.3F);solver.set_interpolated_wall_link(x,8,4,0.7F);}}
        solver.step(17);std::vector<double> out;
        for(std::size_t i=0;i<120;++i){out.push_back(solver.density()[i]);out.push_back(solver.velocity_x()[i]);out.push_back(solver.velocity_y()[i]);}
        emit("lbm_"+std::to_string(scenario),out);
    }
    for(int scenario=0;scenario<9;++scenario){
        cfd::fdtd::Maxwell1DConfig cfg;cfg.cells=96;cfg.courant=0.7;cfg.pml_cells=10;
        if(scenario==1)cfg.boundary=cfd::fdtd::Boundary1D::mur1;
        if(scenario==2)cfg.boundary=cfd::fdtd::Boundary1D::pml;
        if(scenario==3||scenario==7)cfg.boundary=cfd::fdtd::Boundary1D::cpml;
        cfd::fdtd::Maxwell1D solver(cfg);
        if(scenario==4)solver.set_debye_material(38,62,2.1,0.4,2e-11,0.001);
        if(scenario==5)solver.set_drude_material(38,62,2.1,1.4e10,4e8,0.001);
        if(scenario==6)solver.set_lorentz_material(38,62,2.1,0.4,1.4e10,4e8,0.001);
        if(scenario==7)solver.set_tfsf_source(22,[](double t){const double x=(t-8e-11)/2e-11;return std::exp(-0.5*x*x);});
        else solver.initialize_gaussian(0.3,0.045);
        if(scenario==8)solver.set_parallel_lumped_rlc(45,1e-3,1e-6,100.0,1e-9,1e-13);
        solver.step(105);std::vector<double> out(solver.electric().begin(),solver.electric().end());
        out.insert(out.end(),solver.magnetic().begin(),solver.magnetic().end());out.insert(out.end(),solver.polarization().begin(),solver.polarization().end());out.push_back(solver.energy());
        emit("maxwell_"+std::to_string(scenario),out);
    }
    for(int scenario=0;scenario<4;++scenario){
        cfd::fvm::Compressible1DConfig cfg;cfg.cells=40;cfg.periodic=scenario>=2;
        if(scenario%2)cfg.reconstruction=cfd::fvm::CompressibleReconstruction::characteristic_weno5;
        if(scenario%2)cfg.time_integrator=cfd::fvm::CompressibleTimeIntegrator::ssprk3;
        cfd::fvm::CompressibleEuler1D solver(cfg);
        if(scenario<2)solver.initialize({1.0,0.0,1.0},{0.125,0.0,0.1});
        else{std::vector<cfd::fvm::EulerPrimitive1D> q(40);for(std::size_t i=0;i<40;++i)q[i]={1.0+0.1*std::sin(6.283185307179586*(static_cast<double>(i)+0.5)/40.0),0.5,1.0};solver.initialize_profile(q);}
        for(int i=0;i<9;++i)solver.step(0.001);
        std::vector<double> out;for(const auto& u:solver.conserved()){out.push_back(u.density);out.push_back(u.momentum);out.push_back(u.energy);}emit("euler_"+std::to_string(scenario),out);
    }
    cfd::fvm::Diffusion2D diffusion({9,8,1.5,0.7,0.3});diffusion.set_source(1.3);diffusion.set_dirichlet(1.0,2.0,0.0,0.5);(void)diffusion.iterate(45);emit("diffusion2d",{diffusion.field().begin(),diffusion.field().end()});
    {
        cfd::fem::Poisson2DConfig cfg;cfg.diffusivity=1.7;
        cfd::fem::Poisson2D p(cfd::fem::make_rectangle_tri_mesh(5,4,1.3,0.8),cfg);
        p.solve([](cfd::fem::Node2 x){return 2.0+x.x*x.y;},[](cfd::fem::Node2 x){return 0.4+x.x-0.2*x.y;});emit("poisson2d",p.solution());
        cfd::fem::Poisson3D p3(cfd::fem::make_box_tet_mesh(3,3,3));p3.solve([](cfd::fem::Point3 x){return 1.0+x.x;},[](cfd::fem::Point3 x){return 0.2+x.x-x.y;});emit("poisson3d",p3.solution());
        cfd::fem::Heat2DConfig hc;hc.conductivity=0.8;hc.volumetric_heat_capacity=4.0;hc.dt=0.002;
        cfd::fem::Heat2D h(cfd::fem::make_rectangle_tri_mesh(5,5),hc);h.initialize([](cfd::fem::Node2 p){return std::sin(3.141592653589793*p.x)*std::sin(3.141592653589793*p.y);});h.run(4,[](cfd::fem::Node2,double){return 0.3;});emit("heat2d",h.temperature());
    }
    for(int mode=0;mode<2;++mode){
        auto mesh=cfd::fem::make_rectangle_tri_mesh(5,4);cfd::fem::Elasticity2DConfig c;c.youngs_modulus=3.0;c.mode=mode?cfd::fem::ElasticityMode2D::planeStrain:cfd::fem::ElasticityMode2D::planeStress;cfd::fem::Elasticity2D s(mesh,c);
        s.set_dirichlet([](cfd::fem::Node2 p){return p.x==0.0;},[](cfd::fem::Node2){return cfd::fem::Displacement2{};});std::vector<double> t;for(auto p:mesh.nodes)t.push_back(310.0+5.0*p.x);s.solve_thermal(t,0.002,300.0,[](cfd::fem::Node2){return cfd::fem::Displacement2{0.1,-0.2};});std::vector<double> out;for(auto d:s.displacement()){out.push_back(d.x);out.push_back(d.y);}emit("elasticity_"+std::to_string(mode),out);
    }
    {
        cfd::fem::BarModal1DConfig c;c.elements=16;c.young_modulus=100.0;c.density=4.0;cfd::fem::BarModal1D b(c);std::vector<double> out;for(auto m:b.solve(2)){out.push_back(m.angular_frequency);out.insert(out.end(),m.displacement.begin(),m.displacement.end());}emit("modal_bar",out);
        cfd::fem::PiezoelectricBar1D p;p.solve_receive(1.2);out=p.displacement_m();out.insert(out.end(),p.potential_v().begin(),p.potential_v().end());out.push_back(p.electrode_charge_c());p.solve_transmit(10.0);out.insert(out.end(),p.displacement_m().begin(),p.displacement_m().end());out.insert(out.end(),p.potential_v().begin(),p.potential_v().end());out.push_back(p.electrode_charge_c());emit("piezo_bar",out);
        cfd::fem::BilinearPlasticity1D law(200.0,2.0,10.0);out.clear();for(double strain:{0.005,0.03,0.025,-0.03,0.0}){auto r=law.update(strain);out.insert(out.end(),{r.stress,r.tangent,r.plastic_strain,r.equivalent_plastic_strain});}auto r=cfd::fem::compressible_neo_hookean({1.1,0.95,1.03},7.0,11.0);out.push_back(r.energy_density);out.insert(out.end(),r.nominal_stress.begin(),r.nominal_stress.end());emit("material_laws",out);
    }
    {
        std::vector<double> out;auto append=[&](const cfd::rf::Matrix2C& m){for(auto z:{m.a11,m.a12,m.a21,m.a22}){out.push_back(z.real());out.push_back(z.imag());}};
        const cfd::rf::TransmissionLine line{63.0,{0.07,15.2},0.08};const auto s=line.s_parameters(50.0);append(s);append(cfd::rf::s_to_z(s));append(cfd::rf::s_to_abcd(s));append(cfd::rf::cascade_s(s,s));
        auto m=cfd::rf::microstrip_quasi_static(0.001,0.0008,4.2,1.7e9);out.insert(out.end(),{m.characteristic_impedance,m.effective_permittivity,m.phase_velocity,m.guided_wavelength});
        auto c=cfd::rf::coax_quasi_static(0.0005,0.0017,2.1,1.7e9);out.insert(out.end(),{c.characteristic_impedance,c.phase_velocity,c.propagation_constant_rad_per_m,c.guided_wavelength});
        auto w=cfd::rf::rectangular_waveguide_te10(0.02286,0.01016,1.0,1.0,10e9);out.insert(out.end(),{w.cutoff_frequency_hz,w.propagation_constant_rad_per_m,w.guide_wavelength,w.wave_impedance});
        auto cp=cfd::rf::coplanar_waveguide_quasi_static(0.001,0.0002,4.2,1.7e9);out.insert(out.end(),{cp.characteristic_impedance,cp.effective_permittivity,cp.phase_velocity,cp.guided_wavelength});emit("rf_network",out);
    }
    {
        std::vector<double> out;cfd::optics::GaussianBeam b(532e-9,0.0007,0.012,1.0);b.propagate(0.03);b.thin_lens(0.08);out.insert(out.end(),{b.q().real(),b.q().imag(),b.radius(),b.curvature_radius()});
        cfd::optics::MaterialCatalog materials;out.push_back(materials.at("N-BK7").refractive_index_nm(532.0));
        auto f=cfd::optics::fresnel_dielectric(1.0,1.52,0.7);out.insert(out.end(),{f.cos_transmitted,f.rs,f.rp,f.ts,f.tp,f.reflectance_s,f.reflectance_p});
        out.push_back(cfd::optics::thin_film_reflectance_normal(1.0,1.52,532.0,{{1.38,93.0},{2.1,63.0}}));
        cfd::optics::SequentialOpticalSystem system;system.add_surface({cfd::optics::SurfaceType::sphere,0.01,0.05,0.005,1.5});
        for(double x:{-0.002,0.0,0.002}){auto r=system.trace({{x,0.0,0.0},{0.0,0.0,1.0},532.0});out.insert(out.end(),{r.ray.origin.x,r.ray.origin.z,r.ray.direction.x,r.ray.direction.z,r.valid?1.0:0.0});}emit("optics",out);
    }
    {
        using namespace cfd::chemistry;ReactionNetwork network({{"A",0.01,0},{"B",0.01,0}});ElementaryReaction r;r.reactants={{0,1.0}};r.products={{1,1.0}};r.forward={12.0,0.0,0.0};network.add_reaction(r);
        std::vector<double> out=network.source_terms(std::vector<double>{2.0,0.0},300.0);std::vector<double> concentration{1.0,0.0};auto result=integrate_isothermal(network,concentration,300.0,0.2);out.insert(out.end(),concentration.begin(),concentration.end());out.push_back(result.time_advanced);
        IdealGasPhase phase({{"A",10.0},{"B",20.0}},1e5);auto mu=phase.chemical_potentials(std::vector<double>{0.3,0.7},310.0,2e5);out.insert(out.end(),mu.begin(),mu.end());auto p=equilibrate_binary_salt(3.0,2.0,0.0,0.25);out.insert(out.end(),{p.dissolved_a,p.dissolved_b,p.solid_amount,p.reaction_extent,p.saturation_ratio});emit("chemistry",out);
    }
    {
        cfd::acoustics::KSpace2DConfig cfg;cfg.nx=16;cfg.ny=16;cfg.pml_cells=2;
        cfd::acoustics::KSpaceAcoustic2D s(cfg);s.set_uniform_medium(1200.0,998.0);s.set_uniform_nonlinearity(3.0);s.set_power_law_absorption({0.3,1.5,1e6});s.add_sensor(136);s.add_sensor(140);s.initialize_gaussian_pressure(0.5,0.5,0.0002,700.0);s.add_pressure_source(100,[](double t){return 2.0*std::sin(1e6*t);});s.step(4);
        std::vector<double> out=s.pressure();out.insert(out.end(),s.velocity_x().begin(),s.velocity_x().end());out.insert(out.end(),s.velocity_y().begin(),s.velocity_y().end());for(const auto& sensor:s.sensors())out.insert(out.end(),sensor.pressure.begin(),sensor.pressure.end());emit("acoustics",out);
        std::vector<std::size_t> cells{136,140};std::vector<std::vector<double>> traces{s.sensors()[0].pressure,s.sensors()[1].pressure};emit("acoustic_time_reversal",s.time_reversal_reconstruct(cells,traces,true));
    }
    std::cout<<"}\n";
}
