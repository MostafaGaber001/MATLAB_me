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
AF=108 %mm  % outside diameter of flange 
Sa=448 %MPa % ASME SEC2 P.561 row 17 SA-449
Sb=316 %MPa %
Sf=164 %MPa
N_bolts=4
go=5
h=2.5*go
g1=h/3+go
Gi=34+5*2 %assumed  inner 
Go=34+5*2+2*g1 %assumed  and outer diamters of gasket 
C=80  %assumed bolt circle diameter
B=34+5*2  %outside shell diameter or inner 
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

%integral flange stresses laws (INTEGRAL ONLY)
d=(U/V)*ho*go^2
ho=(B*go)^0.5
e=F/ho
L=((t_flange*e+1)/T)+t_flange^3/d
SH=(f*Mo)/(L*g1^2*B)
SR=((1.33*t_flange*e+1)*Mo)/(L*t_flange^2*B)
ST=(Y*Mo)/(t_flange^2*B)-Z*SR
%{d=(U/VL)*ho*go^2
SH_EQ1=SH==Sf
SR_EQ1=SR==Sf
ST_EQ1=ST==Sf
t_SH=double(solve(SH_EQ1,t_flange))
t_SR=double(solve(SR_EQ1,t_flange))
t_ST=double(solve(ST_EQ1,t_flange))