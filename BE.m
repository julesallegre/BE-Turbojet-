clear 
close all


%donnÃ©es Cruise condition

%flight conditions
l=10.700;%Km
M=0.8;
Cp_a=1004;
gamma_a=1.4;
r_a=(Cp_a*(gamma_a-1))/gamma_a;%J/Kg-K

%Componente characteristics:
Fan_pr=1.75;
TPR_IPC=1.3;
TPR_HPC=12;
T_t4=1540;
Eff_comb=0.995;
h_25=-43124e3;%J/kg
Cpg=1147;
gamma_g=4/3;
r_g=(Cpg*(gamma_g-1))/gamma_g;%J/Kg-K

%Losses 
 Inlet_loss = 0.96;
 Fan_is_eff = 0.91;
 IP_comp_is_eff = 0.815;
 Lo_comp = 0.99;
 HPC_is_eff= 0.89;
 Tot_p_loss_cc = 0.97;
 Poly_eff_HPt= 0.85;
 IntT_loss = 0.98;
 Poly_eff_LPt= 0.85;
 Mec_eff_HPshaft = 0.99;
 Mec_eff_LPshaft = 0.99;

%Mass flow and bleed rates :
 mred = 400;% REDUCED MASSFLOW in SECTION 2 Kg/s
 BPR = 5.5;
 Mecpower_HPshaft = 50e3;% W
 Cool_m_HPdist = 0.06; % of overall massflow in section 25, withdrawn at the exit of HP compressor
 Cool_leak_HPT = 0.04; % of overall massflow in section 25, withdrawn at the exit of HP compressor
 Cool_leak_LPT=0;


%CRUISE


%Tables creation
Section={'atm';'1';'2';'13';'16';'18';'21';'24';'25';'3';'31';'4';'41';'44';'45';'5';'6';'8'};
p0=zeros(18,1);
t0=zeros(18,1);
Massflow=zeros(18,1);
tom = table(p0,t0,Massflow);
tom.Properties.RowNames=Section;
tom.Properties.VariableUnits={'Pa','K','Kg/s'};
in=zeros(11,1);
others=table(in);% If the question is not a value --> 1==yes; 0==no; 
others.Properties.RowNames={'Thrust Primary Nozzle';'Primary Nozzle Chocked?';'Static Pressure Primary Nozzle';'Primary Nozzle Area';'Thrust Secondary Nozzle';'Secondary Nozzle Chocked?';'Static Pressure Secondary Nozzle';'Secondary Nozzle Area';'Thrust';'sfc';'Propulsive Efficiency'};



%atmosphere
tom.t0('atm')=288.15-6.5*l;
tom.p0('atm')=101325*(1-0.0225577*l)^(5.25588);

%Air entrance conduit (Point 1)
tom.p0('1')=tom.p0('atm')*(1+(gamma_a-1)*0.5*M^2)^(gamma_a/(gamma_a-1));%isentropique (??)
tom.t0('1')= tom.t0('atm')*(1+(gamma_a-1)*0.5*M^2);

%Point 2
tom.p0('2')=tom.p0('1')*Inlet_loss;
tom.t0('2')=tom.t0('1');

%Flight velocity
V1=M*sqrt(gamma_a*r_a*tom.t0('2'));
Va=V1;

%True massflow
m=(mred*(tom.p0('2')/101325))/sqrt(tom.t0('2')/288.15);
tom.Massflow('2')=m;

%back to point 1
tom.Massflow('1')=m;

%Secodary Flow

%Point 13
tom.Massflow('13')=m/(1+(1/BPR));
tom.p0('13')=Fan_pr*tom.p0('2');
t013_is=tom.t0('2')*(Fan_pr^((gamma_a-1)/gamma_a));
tom.t0('13')=((t013_is-tom.t0('2'))/ Fan_is_eff)+tom.t0('2');

%Secondary nozzle

%Point 16
tom.p0('16')=tom.p0('13');
tom.t0('16')=tom.t0('13');
tom.Massflow('16')=tom.Massflow('13');

%Point 18
tom.t0('18')=tom.t0('16');
tom.p0('18')=tom.p0('16');
tom.Massflow('18')=tom.Massflow('16');

%Critical  conditions (M=1)
pcr18=tom.p0('16')*(1-((gamma_a-1)/(gamma_a+1)))^(gamma_a/(gamma_a-1));
tcr18=2*tom.t0('16')/(gamma_a+1);

if (pcr18>=tom.p0('atm'))
	others.in('Secondary Nozzle Chocked?')=1;
    others.in('Static Pressure Secondary Nozzle')=pcr18;
    t18=tcr18;%Static Temperature Secondary Nozzle
    v18=sqrt(gamma_a*r_a*t18);%==a18
    p018_test=pcr18*(1+(v18^2/(2*Cp_a*t18)))^(gamma_a/(gamma_a-1));% p0('18')==p0('16'), just to test
    
	%Secondary Nozzle Area
	others.in('Secondary Nozzle Area')=(tom.Massflow('18')*r_a*t18)/(pcr18*v18);

	%Thrust Secondary Nozzle: chocked so addition of the term A(p-patm)
	others.in('Thrust Secondary Nozzle')=(tom.Massflow('18')*v18 - tom.Massflow('13')*Va)+(others.in('Secondary Nozzle Area')*(pcr18-tom.p0('atm')));

else
    others.in('Secondary Nozzle Chocked?')=0;
    others.in('Static Pressure Secondary Nozzle')=tom.p0('atm');%total pressure in atm == static, because no speed!!
   
    %static quatities
    p18=others.in('Static Pressure Secondary Nozzle');
    M18=sqrt((2/(gamma_a-1))*((((tom.p0('18')/p18)^((gamma_a-1)/gamma_a))-1)*(2/(gamma_a-1))));
    t18=tom.t0('18')/(1+((gamma_a-1)/2)*M18^2);
    v18=M18*sqrt(gamma_a*r_a*t18);

    %Primary  Nozzle Area
    others.in('Secondary Nozzle Area')=(tom.Massflow('18')*r_a*t8)/(p18*v18);

    %Thrust Primary  Nozzle: not chocked so term A(p-patm)=0
    others.in('Thrust Secondary Nozzle')=(tom.Massflow('18')*v18 - m*Va);
end

%Primary Flow

%Point 21
tom.Massflow('21')=m-tom.Massflow('13');
tom.p0('21')=Fan_pr*tom.p0('2');
t021_is=tom.t0('2')*(Fan_pr^((gamma_a-1)/gamma_a));
tom.t0('21')=((t013_is-tom.t0('2'))/ Fan_is_eff)+tom.t0('2');

%Point 24
tom.p0('24')=TPR_IPC*tom.p0('21');
t24is=tom.t0('21')*(TPR_IPC^((gamma_a-1)/gamma_a));
tom.t0('24')=tom.t0('21')+((t24is-tom.t0('21'))/IP_comp_is_eff);
tom.Massflow('24')=tom.Massflow('21');

%Point 25 HP compressor inlet
tom.p0('25')=Lo_comp*tom.p0('24');
tom.t0('25')=tom.t0('24'); %no work, no heat, dh0=0
tom.Massflow('25')=tom.Massflow('24');

%Point 3 HP compressor outlet
tom.p0('3')=TPR_HPC*tom.p0('25');
t03is=tom.t0('25')*(TPR_HPC^((gamma_a-1)/gamma_a));
tom.t0('3')=tom.t0('25')+((t03is-tom.t0('25'))/HPC_is_eff);
tom.Massflow('3')=tom.Massflow('25');

%rq: use gamma_g and Cpg
%Point 31 Combustion Chamber inlet
tom.p0('31')=tom.p0('3');
tom.t0('31')=tom.t0('3');
tom.Massflow('31')=tom.Massflow('3')-((Cool_m_HPdist+Cool_leak_HPT+Cool_leak_LPT)*tom.Massflow('25'));

%Point 4 Combustion Chamber outlet 
tom.p0('4')=Tot_p_loss_cc*tom.p0('31');%Pressure loss in the CC
tom.t0('4')=T_t4;

%addition of fuel
fth=((Cp_a*(tom.t0('31')-298.15))+(Cpg*(298.15-tom.t0('4'))))/(h_25+(Cpg*(tom.t0('4')-298.15))); %theoric necessary quantity of fuel
f=fth/Eff_comb; %actual quantity of fuel/ quantity of air
tom.Massflow('4')=tom.Massflow('31')*(1+f);

% %Point 41 distributor outlet/rotor inlet - first stage of the turbine
 % stator --> no work, no loss --> total pressure conserved
 tom.p0('41')=tom.p0('4');
 %introduction of cool massflow
 mcool=Cool_m_HPdist*tom.Massflow('25');
 tom.Massflow('41')=tom.Massflow('4')+mcool;
 %Heat diffusion by the cool massflow
 tom.t0('41')=(Cpg*tom.Massflow('4')*tom.t0('4')+Cp_a*mcool*tom.t0('25'))/(Cpg*tom.Massflow('4')+Cp_a*mcool);


%Point 44 - outlet of HP turbine
 %before restitution of cooling air
 %Power Balance:
 Pc=tom.Massflow('25')*Cp_a*(tom.t0('3')-tom.t0('25'));
 Pt=-(Pc/Mec_eff_HPshaft)-Mecpower_HPshaft;
 T044_b=(Pt/(tom.Massflow('41')*Cpg))+tom.t0('41');
 %Polytropic tranformation
 tom.p0('44')=tom.p0('41')*((T044_b/tom.t0('41'))^(gamma_g/(Poly_eff_HPt*(gamma_g-1))));
 %after restitution of the cooling air
 mcoolhpt=Cool_leak_HPT*tom.Massflow('25');
 tom.Massflow('44')=tom.Massflow('41')+mcoolhpt;
 %Heat diffusion
 tom.t0('44')=(Cpg*tom.Massflow('41')*T044_b+Cp_a*mcoolhpt*tom.t0('25'))/(Cpg*tom.Massflow('41')+Cp_a*mcoolhpt);

%Point 45 - inlet Low Pressure Turbine
 tom.p0('45')=IntT_loss*tom.p0('44');
 %no work, no heat
 tom.t0('45')=tom.t0('44');
 tom.Massflow('45')=tom.Massflow('44');

%Point 5 - LP Turbine outlet
 %before restitution of cooling air
 %Power balance (no withdrawn power)
 plpc=tom.Massflow('21')*Cp_a*(tom.t0('24')-tom.t0('21'));
 plpt=-plpc;
 T05_b=(plpt/(tom.Massflow('45')*Cpg))+tom.t0('45');
 %Polytropic transformation
 tom.p0('5')=tom.p0('45')*((T05_b/tom.t0('45'))^(gamma_g/(Poly_eff_LPt*(gamma_g-1))));
 %after restitution of the cooling air(A17) 
 mcoollpt=Cool_leak_LPT*tom.Massflow('25');
 tom.Massflow('5')=tom.Massflow('45')+mcoollpt;
 %Heat diffusion
 tom.t0('5')=(Cpg*tom.Massflow('45')*T05_b+Cp_a*mcoollpt*tom.t0('25'))/(Cpg*tom.Massflow('45')+Cp_a*mcoollpt);

%Primary Nozzle

%Point 6 - Primary nozzle inlet
 %No work, no heat and no losses in the nozzle (simple adiabatic conduit)
 tom.p0('6')=tom.p0('5');
 tom.t0('6')=tom.t0('5');
 tom.Massflow('6')=tom.Massflow('5');

%Point 8 -  Primary Nozzle outlet
 tom.t0('8')=tom.t0('6'); %no work, no heat
 tom.p0('8')=tom.p0('6');
 tom.Massflow('8')=tom.Massflow('6');

 %Critical  conditions (M=1)
 pcr8=tom.p0('6')*(1-((gamma_g-1)/(gamma_g+1)))^(gamma_g/(gamma_g-1));
 tcr8=2*tom.t0('6')/(gamma_g+1);

 if (pcr8>=tom.p0('atm'))
    others.in('Primary Nozzle Chocked?')=1;
    others.in('Static Pressure Primary Nozzle')=pcr8;
    t8=tcr8;%Static Temperature Primary Nozzle
    v8=sqrt(gamma_g*r_g*t8);%==a18
    p08_test=pcr8*(1+(v8^2/(2*Cpg*t8)))^(gamma_g/(gamma_g-1));% p0('8')==p0('6'), just to test
    
    %Primary  Nozzle Area
    others.in('Primary Nozzle Area')=(tom.Massflow('8')*r_g*t8)/(pcr8*v8);

    %Thrust Primary  Nozzle: chocked so addition of the term A(p-patm)
    others.in('Thrust Primary Nozzle')=(tom.Massflow('8')*v8 - tom.Massflow('21')*Va)+(others.in('Primary Nozzle Area')*(pcr8-tom.p0('atm')));

 else
    others.in('Primary Nozzle Chocked?')=0;
    others.in('Static Pressure Primary Nozzle')=tom.p0('atm');%total pressure in atm == static, because no speed!!
    
    %static quatities
    p8=others.in('Static Pressure Primary Nozzle');
    M8=sqrt((2/(gamma_g-1))*((((tom.p0('8')/p8)^((gamma_g-1)/gamma_g))-1)*(2/(gamma_g-1))));
    t8=tom.t0('8')/(1+((gamma_g-1)/2)*M8^2);
    v8=M8*sqrt(gamma_g*r_g*t8);

    %Primary  Nozzle Area
    others.in('Primary Nozzle Area')=(tom.Massflow('8')*r_g*t8)/(p8*v8);

    %Thrust Primary  Nozzle: not chocked so term A(p-patm)=0
    others.in('Thrust Primary Nozzle')=(tom.Massflow('8')*v8 - m*Va);
end

%Global Characteristics

%Total Thrust
 others.in('Thrust')=others.in('Thrust Primary Nozzle')+others.in('Thrust Secondary Nozzle');

%Propulsive efficiency
 others.in('Propulsive Efficiency')=2*Va*((tom.Massflow('8')*v8)+(tom.Massflow('18')*v18)-(tom.Massflow('1')*Va))/((tom.Massflow('8')*v8^2)+(tom.Massflow('18')*v18^2)-(tom.Massflow('1')*Va^2));

%Specific Consumption
others.in('sfc')=(tom.Massflow('31')*f)/others.in('Thrust');


tom
others









