%{
Author:Mostafa_Gaber
Date_Created:5/5/2024
============================================================
=========Pressure vessel calculations based on div1=========
============================================================
%}
clc 
clear
%{
============================================================
============Pressure vessel size calculation================
============================================================
%}
%===================Reactants parameters====================
Oil_liters=22          %Litre    % 15L of oil to be used   
Alcohol_to_Oil=28/1              % Methanol to oil molar ratio = 28/1
Alcohol_molar_vol=4.05*10^(-5)   % 1mol --> 4.05*10^(-5)m^3
                       %m3/mol
Alcohol_density=792    %kg/m^3   % Methanol density 792 kg/m^3 
Oil_density=903.1      %kg/m3    % 0.9031 g/mL--> 903.1 kg/m3 waste canola oil density
Oil_molar_mass=876.6   %kg/kmol  % 876.6 kg/kmol waste canola oil
%Cat_density=12020     %kg/m3    % pd/al2o3 12.02 g/cm^3 --> g/mL -> 12020 kg/m3
%Cat_percent=5                   % 5% of oil weight
%=============== Oil,Alcohol and Catlayst vols ============
Oil_vol=Oil_liters*10^(-3)                         %m^3
Oil_mass=Oil_density*Oil_vol                       %kg     
Oil_moles= Oil_mass/Oil_molar_mass                 %kmol    
Alcohol_moles=Alcohol_to_Oil*Oil_moles             %kmol
Alcohol_vol=Alcohol_moles*10^(3)*Alcohol_molar_vol %m^3
%Cat_mass=(Cat_percent/100)*(Oil_mass)              %kg
%Cat_vol=Cat_mass/Cat_density                       %m^3
%=========================Nitrogen vol=======================
N2_density=1.2  %kg/m^3          % Nitrogen density
T1=298          %k°              % T1 is temp. at which N2 will be introduced inside
R=0.0821   %Litre*atm/mol*k°     % gas constant
T2=543          %k°              % T2 temperature at which nitrogen operate
a=1.39     %Litre^2*atm/mol^2    % van der waal const. 
b=0.0391     %Litre/mol          % van der waal const.
P2_operating=10/0.101325  %atm   % P2 is pressure at which N2 will operate
n=linspace(0,50,150) ;    %moles              % number of moles of N2
V=linspace(0,50,150);     %Litres             % volume of gas of N2
[V,n]=meshgrid(V,n);        
P2=((n*R*T2)./(V-n*b))-(a*n.^2)./(V.^2) % P2 pressure at which nitrogen operates
[row,col]=find(P2>P2_operating & P2<100)
%P2=P2.*(P2>=P2_required)    
%P2=P2.*(P2<200)
N2_vol=V(row(1),col(1))    % Litres
N2_moles=n(row(1),col(1))  % moles
P1=((N2_moles*R*T1)/(N2_vol))*0.101325 %MPa  % P1 is pressure at which N2 will be introduced inside
PV_Size=Oil_vol+Alcohol_vol+N2_vol*10^(-3)%Cat_vol  %m^3 %final size of PV
%{
============================================================
==========Thickness of cylindrical shell and Head=========== 
=========under internal pressure and Content weight=========
============================================================
%}
%======================PV Parameters========================
PV_Material_Yield=142 %MPa   % Stainless steel SA–790 yield strength
                             % from ASME SEC.2 Table 1A p.86 @ 270 °C row 24
%PV_Size=50*10^(-3)   %m^3   % PV size ceiled from 47.7L to 50L
PV_Length=80*10^(-2)  %m     % PV Length assumed 60 cm
syms PV_In_Dia               % PV inner diameter
PV_Elipsoidal_vol= (pi/24)*PV_In_Dia^3 % ASME Elipsoidal head volume formula , Elipsoidal is a 2:1 Eliptical head
eqn=PV_In_Dia==((PV_Size-PV_Elipsoidal_vol)*4/(PV_Length*pi))^0.5
PV_In_Dia=real(double(solve(eqn,PV_In_Dia))) %m %inner diameter of the PV
PV_Elipsoidal_vol=double(subs(PV_Elipsoidal_vol)) %m^3 
%====================Stresses calculations==================
content_mass=Oil_mass+Alcohol_vol*Alcohol_density%Cat_mass %kg % N2 is neglected
content_force=content_mass*9.81
%======Assuming thin-walled=================================
Shell_thickness=[]
FS=2                                    % Factor of safety 
PV_In_Dia=PV_In_Dia*10^3           %mm  % turning it to mm
P2_operating=P2_operating*0.101325 %MPa % turning it to MPa
syms t_thin
Sigma_hoop=(P2_operating*PV_In_Dia/2)/t_thin
Sigma_axial=(P2_operating*PV_In_Dia/2)/(2*t_thin)+content_force/((pi/4)*(PV_In_Dia+2*t_thin)^2-(pi/4)*PV_In_Dia^2)
t_thin_hoop=double(solve(Sigma_hoop==(PV_Material_Yield/FS),t_thin))
t_thin_axial=double(solve(Sigma_axial==(PV_Material_Yield/FS),t_thin))
t_thin_axial=max(t_thin_axial)
if (PV_In_Dia/max(t_thin_hoop,t_thin_axial))>10 %checking thin wall condition
Shell_thickness=max(t_thin_hoop,t_thin_axial)
else
%the condition D/t>10 was not satisfied so the assumption of thin walled
%not valid
%=====Using Div1 Equations=================================
E=1                    % joint efficiency
t_hoop=(P2_operating*PV_In_Dia/2)/((PV_Material_Yield*E/FS)-0.6*P2_operating)
t_axial=(P2_operating*PV_In_Dia/2)/(2*PV_Material_Yield*E/FS+0.4*P2_operating)
Shell_thickness=max(t_hoop,t_axial)
end
%=====Thickness of ASME Elipsoidal head====================
E=1                    % joint efficiency
a_head=PV_In_Dia/2 %mm % From geometry
b_head=a_head /2   %mm % a 2:1 Elipsoidal head i.e. a/b =2
K=(1/6)*(2+(a_head/b_head)^2)
t_elip=(K*P2_operating*PV_In_Dia)/(2*PV_Material_Yield*E/FS-0.2*P2_operating) %mm
%{
============================================================
==========Reinforcement of Shell and Formed head============
============================================================
%}
%=====Parameters Definitions===============================
A = [];  % total cross‐sectional area of reinforcement
         % required 
A1= [];  % area in excess thickness in the vessel wall
         % available for reinforcement 
A2= [];  % area in excess thickness in the nozzle wall
         % available for reinforcement 
A3 =[];  % area available for reinforcement when the
         % nozzle extends inside the vessel wall (see
A5 =[];  % cross‐sectional area of material added as re-
         % inforcement 
A41 =[]; A42=[];A43 =[]; % cross‐sectional area of various welds available for reinforcement 
c =[];          % corrosion allowance
D =[];          % inside shell diameter
Dp =[];         % outside diameter of reinforcing element (ac-
                % tual size of reinforcing element may exceed
                % the limits of reinforcement established by
                % UG-40; however, credit cannot be taken for
                % any material outside these limits)
d =[];          % finished diameter of circular opening
                % or finished dimension

E =[];          % (see definitions for tr and trn)
E1 =[];         % joint efficiency "carefull definitions changed from above"
F =[];          % correction factor 
fr =[];         % strength reduction factor, not greater than
                % 1.0 [see UG-41(a)] "ratio between reinforcement allowable stress and vessel strength"
fr1 =[];        % Sn/Sv for nozzle wall inserted through the vessel wall
                % = 1.0 for nozzle wall abutting the vessel wall
fr2 =[];        % Sn/Sv
fr3 =[];        % (lesser of Sn or Sp)/Sv
fr4 =[];        % Sp/Sv
h =[]           % distance nozzle projects beyond the inner
                % surface of the vessel wall.
K1 =[];         % spherical radius factor (see definition of tr
                % and Table UG-37)
L =[];          % length of projection defining the thickened
                % portion of integral reinforcement of a nozzle
                % neck beyond the outside surface of the
P =[];          % internal design pressure (see UG-21), psi (MPa)
R =[];          % inside radius of the shell course under consideration
Rn =[];         % inside radius of the nozzle 
S =[];          % allowable stress value in tension 
Sn =[];         % allowable stress in nozzle, psi (MPa) (see S above)
Sp =[];         % allowable stress in reinforcing element (plate), psi (MPa) (see S above)
Sv =[];         % allowable stress in vessel, psi (MPa) (see S above)   
t =[];          % specified vessel wall thickness which thickness needed due to stress
                % +corroision allowance + tolarence
te =[];         % thickness or height of reinforcing element
ti =[];         % nominal thickness of internal projection of nozzle wall
tn =[];         % nozzle wall thickness.29 Except for pipe, this
                % is the wall thickness not including forming
tr =[];         % required thickness of a seamless shell based
                % on the circumferential stress, or of a formed
                % head, computed by the rules of this Division
                % for the designated pressure, using E = 1
trn =[];        % required thickness of a seamless nozzle wall,
W =[];          % total load to be carried by attachment welds (see UG-41)
W11=[];         % Weld load on strength path 1-1 look p.52 div1
W22=[];         % Weld load on strength path 2-2 look p.52 div1
weld_leg11=[];  % Weld leg on strength path 1-1 look p.52 div1
weld_leg22=[];  % Weld leg on strength path 2-2 look p.52 div1
weld_leg=[];    % weld leg on case no reinforcement pad
%=====Strengthes===============
P=P2_operating % 10 MPa
Sv=PV_Material_Yield
Sn=100.88 %MPa % nozzle assumed to be forged 
               % Forgings from ASME SEC.2 Table 1A p.86 @270 °C row 25
Sp=100.88 %MPa % Plate from ASME SEC.2 table 1A p.86 @270 °C row 34
Sw=min([Sp,Sn,Sv]) %MPa % Weld allowable strength = minimum of two attached reinforcement
weld_shear_strength=0.75*Sw
%=====Dimensions===============
D=PV_In_Dia
d=34% based on  valve diameter "Need to choose globe valve"       
Rn=d/2
c=2             % corroision
E=1             % assumed factor
E1=1            % assumed factor
%te=0           % we dont have pad reinforcing for now
ti=0            % we have abutting nozzle
trn=(P*Rn)/(Sn*E1-0.6*P) % nozzle required thickness due to stresses
tn=ceil(2*trn)           % nozzle real thickness used 
tr=Shell_thickness      
t=Shell_thickness+3                     %actual thickness used
Dp=2*max(d,(Rn+tn+t))
%====Factors===================D
F=1
fr=1
fr1=1
fr2=Sn/Sv
fr3=min((Sn/Sv),(Sp/Sv))
fr4=Sp/Sv
%=====Calculations=============
A=d*tr*F+2*tn*F*(1-fr1) 
A11=d*(E1*t-F*tr)-2*tn*(E1*t-F*tr)*(1-fr1)
A12=2*(t+tn)*(E1*t-F*tr)-2*tn*(E1*t-F*tr)*(1-fr1)
A1=max(A11,A12)
A21=5*(tn-trn)*fr2*t
A22=5*(tn-trn)*fr2*tn
A2=min(A12,A22)
% A31=5*t*ti*fr2
% A32=5*ti*ti*fr2
% A33=2*h*ti*fr2
% A3=min(A31,A32,A33)   % no area available inside 
W=(A-A1)*Sv       %Newton  % total weld load assumed to be tensile      
weld_leg=(W*1.414)/(weld_shear_strength*pi*d)
A41=(weld_leg)^2*fr2    % Area available in outward weld
% A43=(weld_leg)^2*fr2  % Area available in inward weld
Total_Area_avail=A1+A2+A41%+A43+A3
if Total_Area_avail>A 
    disp('No Reinforcement need')
else  disp('Reinforcement needed')
    syms te 
   %  A21=5*(tn-trn)*fr2*t
     A22=2*(tn-trn)*(2.5*tn+te)*fr2
     A2=A22
     weld_leg11=weld_leg
     weld_leg22=weld_leg
     A41=(weld_leg11)^2*fr3
     A42=(weld_leg22)^2*fr4
     A5=(Dp-d-2*tn)*te*fr4
     eqn1=A5+A2+A41+A42+A1==A
    te= double(solve(eqn1,te))
end
%{
%{
============================================================
================Supports calculation========================
============================================================
%}
%=====Supports Parameters=======
FS_supp=2
no_supports=4
Dist_from_Center=122        %mm
Incline_angle_on_vertical=-2 %degree
Support_height=300          %mm
Support_length=Support_height*cosd(Incline_angle_on_vertical)
Support_material=200          %MPa % yield assumed ====bl7ob=== for now
syms supp_dia
%=====Calculations==============
supp_thick=0.05:0.05:4       %mm       
size(supp_thick,2)           %no of elements in thickness
D_tensile=zeros(1,size(supp_thick,2)) 
D_comp=zeros(1,size(supp_thick,2))
supp_area=((pi/4)*((supp_dia+supp_thick).^2-(supp_dia)^2)) %mm^2 %support area
supp_inertia=(pi/64)*((supp_dia+supp_thick).^4-(supp_dia)^4)%mm^4 %support inertia
Elipsoidal_head_volume=((pi/24)*(PV_In_Dia+2*t)^3)*10^(-9)-PV_Elipsoidal_vol
PV_density=7805 %kg/m3 %SA_790 density
PV_volume=(pi/4)*((PV_In_Dia+2*t)^2-PV_In_Dia^2)*10^(-6)*PV_Length+Elipsoidal_head_volume%+cover_volume
% Total_mass=content_mass+(PV_volume)*PV_density %kg
Total_mass=content_mass+(41879059.88*10^-9)*PV_density %kg
Total_force=Total_mass*9.81/no_supports
Moment_At_Ground=Total_force*sind(Incline_angle_on_vertical)*Support_length
Sigma_bending=Moment_At_Ground*((supp_dia+supp_thick)/2)./supp_inertia
Max_compressive=Sigma_bending+Total_force*cosd(Incline_angle_on_vertical)./supp_area
Max_tensile=Sigma_bending-Total_force*cosd(Incline_angle_on_vertical)./supp_area
eqn2=Max_compressive==Support_material/FS_supp
eqn3=Max_tensile==Support_material/FS_supp
for i= 1:1:size(supp_thick,2)
D_comp(1,i)=max(double(solve(eqn2(1,i),supp_dia)))    %mm
D_tensile(1,i)=max(double(solve(eqn3(1,i),supp_dia))) %mm
end
%}

%}
%{
============================================================
==============FLANGE design using append.14==============
============================================================
%}
%==========Appednix 14 Parameters defintion============
%Except as given below, the symbols used in the equa-
%tions of this Appendix are defined in 2-3.
AF =[] %outside diameter of flat head and shell
Bn =[] % diameter of central opening (for nozzle, this is in-
       % side diameter and for opening without nozzle,
       % diameter of opening)
Bs =[] % inside diameter of shell (measured below ta-
       % pered hub, if one exists)
E_theta=[] % slope of head with central opening or nozzle
       % times the modulus of elasticity, disregarding
       % the interaction of the integral shell at the outside
       % diameter of the head, psi (MPa)
MH =[] % moment acting at shell‐to‐flat head juncture
P =[]  % internal design pressure (see UG-21)
t =[]  % flat head nominal thickness
%========Appendix 2 parameters==========================
AF =[] %outside diameter of FLANGE
Am=[]% max of Am1 , Am2 
Ab=[]% Area of cross-section area of non threaded area in bolt
G=[] % diameter at location of gasket load reaction i.e mean gasket diameter
Gi=[]% inner diameter of gasket
Go=[]% outer diameter of gasket
m=[] % gasket factor, obtain from Table 2-5.1 [see Note in 2-5(c)(1)]
y =[]  % gasket or joint‐contact‐surface unit seating load, [see Note 1, 2-5(c)]
N= [] % contact length on gasket width
b=[]  % Effective gasket or joint‐contact‐surface seating width [see Note in 2-5(c)(1)]
g1 =[] % thickness of hub at back of flange
go =[] % thickness of hub at small end
Wm1=[] % minimum required bolt load for the operating conditions
Wm2 =[]% minimum required bolt load for gasket seating [see 2-5(c)

B1 =[] % B + g1 for loose type flanges and for integral
       %type flanges that have calculated values h / ho
       %and g1 / go which would indicate an f value of
       %less than 1.0, although the minimum value of f
       %permitted is 1.0.
       %= B + go for integral type flanges when f is equal
       %to or greater than one

Bs =[]% bolt spacing. The bolt spacing may be taken as
      % the bolt circle circumference divided by the
      % number of bolts or as the chord length between
      % adjacent bolt locations.
Bsc =[]% bolt spacing factor
Bsmax =[] % maximum bolt spacing
B =[]% inside diameter of flange. When B is less than
     % 20g1, it will be optional for the designer to sub-
     % stitute B1 for B in the formula for longitudinal
     % stress SH
C =[] % bolt‐circle diameter
c =[] % basic dimension used for the minimum sizing of
      % welds equal to tn or tx, whichever is less
Cb =[]% conversion factor = 0.5 for U.S. Customary calculations; 2.5 for SI
      % calculations
       
       
  F =[]%factor for integral type flanges (from Figure2-7.2)

  Sa=[] %allowable bolt stress at atmospheric temperature (see UG-23)

  Sb =[ ]% allowable bolt stress at design temperature (see UG-23)
Sf =[]   % allowable design stress for material of flange at
         % design temperature (operating condition) or at-
         % mospheric temperature (gasket seating), as
         % may apply (see UG-23)
SH =[] %calculated longitudinal stress in hub       
SR =[] %calculated radial stress in flange   
ST =[] %calculated tangential stress in flange       
Mo =[]% total moment acting upon the flange, for the op-
      % erating conditions or gasket seating as may apply (see 12-4)
V =[] % factor for integral type flanges (from Figure%2-7.3) 
 f =[] % hub stress correction factor for integral flanges
       % from Figure 2-7.6 (When greater than one, this
       % is the ratio of the stress in the small end of hub
       % to the stress in the large end.) (For values below
       % limit of figure, use f = 1.)
ho =[] % =factor=(B*go)^0.5       
L=[]   % factor
e=[]   % factor
K =[] % ratio of outside diameter of flange to inside dia-
      % meter of flange A/B --> AF/B
h =[] % hub length
N_bolts=[]% Number of bolts
%===========Calculations==========================
syms t_flange
t_flange_output=zeros(1,3)
P=10 %MPa
m=3
y=69
% t_flange=28
AF=480 %mm  % outside diameter of flange 
Sa=448 %MPa % ASME SEC2 P.561 row 17 SA-449
Sb=316 %MPa %
Sf=164 %MPa
N_bolts=5
go=9
h=2.5*go
g1=h/3+go
Gi=314 %assumed  inner 
Go=347 %assumed  and outer diamters of gasket 
C=413.5  %assumed bolt circle diameter
B=314  %outside shell diameter or inner 
ho=(B*go)^0.5
G=(Go+Gi)/2
N=(Go-Gi)/2
bo=N/2
b=2.5*(bo)^0.5
Wm1=0.785*(G^2)*P+(2*b*pi*G*m*P)
Wm2=pi*G*y*b
Am2 = Wm2/Sa
Am1 = Wm1/Sb
Am=max(Am1,Am2)
Ab=Am    
W_Operating=Wm1
W_Seating=(Am+Ab)*Sa/2
Mo=W_Seating*(C-G)/2
K=AF/B
[T,U,Y,Z] = Values_of_TUYZ(K)
[F,V,f,FL,VL,fL] = Values_of_FVf(g1,go,h,ho)
%{
%integral flange stresses laws (INTEGRAL ONLY)
d=(U/V)*ho*go^2
ho=(B*go)^0.5
e=F/ho
L=((t*e+1)/T)+t^3/d
SH=(f*Mo)/(L*g1^2*B)

%}
d=(U/VL)*ho*go^2
e=FL/ho
L=((t_flange*e+1)/T)+t_flange^3/d
SH=(fL*Mo)/(L*g1^2*B)
SR=((1.33*t_flange*e+1)*Mo)/(L*t_flange^2*B)
ST=(Y*Mo)/(t_flange^2*B)-Z*SR
SH_EQ1=SH==Sf
SR_EQ1=SR==Sf
ST_EQ1=ST==Sf
t_SH=double(solve(SH_EQ1,t_flange))
t_SR=double(solve(SR_EQ1,t_flange))
t_ST=double(solve(ST_EQ1,t_flange))

%{
============================================================
===================Flat head design=========================
============================================================
%}
C_factor =0.3  % a factor depending upon the method of attachment
       % of head, shell dimensions, and other items as listed
       % in (d) below, dimensionles
       E=1
S_flat=164
hg=(C-G)/2
d_shortspan=G  %mm
%====careful!! value calculated after reinforcing====================
t_flat_operating=d_shortspan*(2*((C_factor*P)/(S_flat*E)+1.9*W_Operating*hg/(S_flat*E*d_shortspan^3)))^0.5
t_flat_Seating=d_shortspan*(2*((C_factor)/(S_flat*E)+(1.9*W_Seating*hg)/(S_flat*E*d_shortspan^3)))^0.5

%{
%{
============================================================
===================Mixer Shaft design=======================
============================================================
%}
%{
============================================================
===================Mixed fluid properties===================
============================================================
%}

Operating_speed=600 %rpm   % SPEED OF AGITATOR 
Mix_density=(Oil_mass+Alcohol_vol*Alcohol_density)/(Oil_vol+Alcohol_vol) %kg/m^3
Oil_viscosity=46.2*0.001 %cP   
Alcohol_viscosity=0.543*10^(-3)
x_alcohol=Alcohol_moles/(Alcohol_moles+Oil_moles)
x_oil=Oil_moles/(Alcohol_moles+Oil_moles)
first_term=(Alcohol_viscosity*x_alcohol*(Alcohol_vol*Alcohol_density/Alcohol_moles)^0.5)
second_term=(Oil_viscosity*x_oil*(Oil_molar_mass)^0.5)
Mix_viscosity=(second_term+first_term)/(x_oil*(Oil_molar_mass)^0.5+x_alcohol*(Alcohol_vol*Alcohol_density/Alcohol_moles)^0.5)



Motor_Power=1000 %watt
Torque=(Motor_Power)/(Operating_speed*(2*pi/60))%N.m
Torque=Torque*10^3
Motor_yield_strength=140
Motor_density=7860 %kg/m3
D_shaft=0:1:45
Motor_CS=(pi/4)*(D_shaft.*10^(-3)).^2
Motor_CS_mm=(pi/4)*(D_shaft).^2
Sigma_t=-10 %MPa
Sigma_r=-10
tao=(Torque*D_shaft/2)./((pi/32)*D_shaft.^4)
Sigma_z=-10
Sigma_von=(3*tao.^2).^0.5
%}


