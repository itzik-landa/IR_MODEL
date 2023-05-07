function DryingSimIR2022
%%
LengthDryer = 2.0;
p0 = 2.4*39e3/LengthDryer;
PostDryer = 0.02;
fig_out = 1;

%%
% mind the v_air -> maibe cancel front to back
% T pressambient
% T ambient
% air dryer
% reduing blanket temperature and hence evaporation rate
% evaporative cooling
% measure the cooling time constant of the blanket
% blanket h graph
% extend the thickness when dry
%%
% Absorbed Energy/Input Energy = 1 - exp( - Alpha*Percent Filling*Depth)
% Alpha=0.276 1/(%*u)
% Percent Filling=7% or 20%
% Depth = x
% X = Humidity_Ratio_PG: has temperature limitation
%%
tic;

%% Dryer Length


Width = 1;% m

% Power_Precentage = 100;%
%%
PercentFilling = 62;%%[7 20 500]

ReleaseThickness = 20;%um[50 5]
ReleaseLayersNumber = 1;%[2]
AbsorbingLayerThickness = 50;%um[150 50]
AbsorbingLayersNumber = 10;%[3 1]%50;
AboveFabricLayerThickness = 200;%[150]
AboveFabricLayersNumber = 1;
FabricLayerThickness = 150;%[100]
FabricLayersNumber = 1;
GripLayerThickness = 100;%[100]
GripLayersNumber = 1;

% p0 = 200e3;%kW/m2[200e3]
% LengthDryer = 1.045;%0.250*8;%

% p0 = 39e3 * Power_Precentage/100 /(0.5*Width)*2.6;%1.45;% [2/0.5 of m/m]3.250
%%
Te = 0.60;
%Rw = 0.15;
%Ew = 1 - Rw;
% A7p5um = 0.065/5 * 7.5;
% Ei = Ee * Ew * A7p5um;
% E_ti = 1 - A7p5um;
Rr = 0.05;
Tr = 1 - Rr;
% Eb_ti = Ee * Ew * E_ti;
Tb_tr = Te * Tr;

% Pi = Ei * P0;
Pi = 0;%
% Pb_ti = Eb_ti * P0;
pb_tr = Tb_tr * p0;


%%


% AreaCoverage_0 = 10;


AreaCoverage_PB(1) = 0;
AreaCoverage_PB(2) = 100;
AreaCoverage_PB(3) = 0;
AreaCoverage_PB(4) = 0;%100;%[0 100]
AreaCoverage_PB(5) = 50;%[50 50]
AreaCoverage_PB(6) = 0;%50;%[100 100]
AreaCoverage_PB(7) = 0;


TotalCoverage = sum(AreaCoverage_PB);

% Test_Time = (1.122+1.878)/Print_Speed;%10; % sec 1:0.2455, 2:

%% Test Matrix
TestNumber = fig_out;%604;
FigureNumber = TestNumber;
NumberOfIPUs = 7;
T_AirDryer = 300;%1200;%2150;%270;%
% ExtraDryerLength = 0.0;%1.80;%1.05;%2.05;%
% LengthDryer = 0.250*8+ExtraDryerLength;
T_IPU_All = 100;%70;%90;%165;%133;%200;%
T_Blanket0 = 68;%70;%80;%
Print_Speed_SPH = 6500;%13000;%10000;%
%%
% PercentFilling=7;%%[7,20]100500
% ReleaseThickness = 50;%um[50]
%%

% ReleaseThickness = 50;%um[50]

dzRelease = ReleaseThickness/ReleaseLayersNumber;
dzAbsorbing = AbsorbingLayerThickness/AbsorbingLayersNumber;
dzAboveFabric = AboveFabricLayerThickness/AboveFabricLayersNumber;
dzFabricLayer = FabricLayerThickness/FabricLayersNumber;
dzGripLayers = GripLayerThickness/GripLayersNumber;
%%
% CalculationLayerThickness_um = 50;%um[50,10]
% LayerThickness = CalculationLayerThickness_um*1e-6;
% dz0 = LayerThickness;%500e-6; %m
% z = (ReleaseThickness+AbsorbingLayerThickness+AboveFabricLayerThickness+...
%     FabricLayerThickness+GripLayerThickness)*1e-6;%500e-6;%1500e-6;%5000e-6; %m
% nz0 = int32(z/dz0); %ul

dz = [dzRelease*ones(ReleaseLayersNumber,1) ; ...
    dzAbsorbing*ones(AbsorbingLayersNumber,1) ; ...
    dzAboveFabric*ones(AboveFabricLayersNumber,1) ; ...
    dzFabricLayer*ones(FabricLayersNumber,1) ; ...
    dzGripLayers*ones(GripLayersNumber,1)]*1e-6;%dz0*ones(nz0,1); %m

Print_Speed = 1.7/6500*Print_Speed_SPH;
PostDryerLength = PostDryer;%0.375;%0;%0.9;%
Length = 0.375*NumberOfIPUs+LengthDryer+PostDryerLength;%-ExtraDryerLength+0.9;%0.375*4+LengthAirDryer+0.1;%Print_Speed * Test_Time; % m
if NumberOfIPUs == 4
    Length = Length + 0.125 + 1;%-ExtraDryerLength;
end


for idx = 1:7
PB_Start(idx) = (0.375*(idx-1)+0.2+0.175/2)/Length;
end

T_IR = 120;%25;% [120];
h_air_IR = 50;%
% q_radiation_Blanket = Pb_ti;%
q_radiation_Blanket = pb_tr;%
Q_radiation_Ink = Pi;%


% TestTitle = '10k 7C 250%:  ';


for idx=1:7
T_IPU (idx) = T_IPU_All;
end

h_IPU = 200;



for idx=1:7
    AirIPU_Start(idx) = (0.375*(idx-1)+0.2)/Length;
    AirIPU_end(idx) = (0.375*(idx-1)+0.2)/Length;
end



h_AirDryer = 200;

AirDryer_Start = 0.375*NumberOfIPUs/Length;
AirDryer_end = (0.375*NumberOfIPUs+LengthDryer)/Length*0;

Length_IR_Dryer = LengthDryer;
Radiation_start = 0.375*NumberOfIPUs/Length;
Radiation_end = (0.375*NumberOfIPUs+Length_IR_Dryer)/Length;
Rad_Ratio = 3/4;
Rad_Number = 1;
Rad_Length = Length_IR_Dryer*Rad_Ratio/Rad_Number;
IR_EffectiveLength = Rad_Ratio;
%%
% AirIPU_Start = (404+1200)/(1122+1878);
% %% 1 2 3 4
% N_IPU = 1;
% AirIPU_end = (404+1200+N_IPU*93-1)/(1122+1878);
% h_IPU = 100;%
% T_IPU = 150;

% LayerNumber = 2:(2+4);%2:(2+0);%8:(8+0);%2:(2+24);%8:(8+24);%
% LayerThickness = 50e-6;%20e-6;%
%%
% AbsorptionLayerStart = ReleaseThickness;%50;%um[50,20]
% AbsorptionLayerEnd = AbsorptionLayerStart + AbsorptionLayerThickness;%um
% NumberOfAbsorptionLayers = AbsorbingLayersNumber;%#
IndexOfFirstAbsorptionLayer = ReleaseLayersNumber + 1;%#
Alpha=0.276;%1/(%*u)[0.276]
for j=1:AbsorbingLayersNumber+1
    AbsorptionDepth = (j-1) * dzAbsorbing;%um
    AbsorbedEnergyDevidedByInputEnergy(j) = 1 - exp(-Alpha*PercentFilling/100*AbsorptionDepth);%%
end
AbsorptionPerLayer = diff(AbsorbedEnergyDevidedByInputEnergy);%%
LayerNumber = IndexOfFirstAbsorptionLayer:(IndexOfFirstAbsorptionLayer+AbsorbingLayersNumber-1);%
q_radiation_Layer = q_radiation_Blanket*AbsorptionPerLayer;
% figure(468);plot(LayerNumber,AbsorptionPerLayer,'.-');
figure(468);plot(linspace(0,AbsorbingLayerThickness,length(AbsorbedEnergyDevidedByInputEnergy)),AbsorbedEnergyDevidedByInputEnergy,'.-');

%%

HelpLayer = 10e-6;

%% 0.2 0.28 0.56

%%

% Test_Time = 1.122/Print_Speed;%10; % sec 1:0.2455, 2:
% Length_IR_Dryer = Length*(Radiation_end-Radiation_start);



%% Other Dryer Setup

%% 0 1.5
v_air = 0;%-Print_Speed;%1.5;%0;%1.5; % m/s

RH_ambient = 50; % %

T_Room_ambient = 25; % %

T_Press_Ambient = 45;%25;% [45];

h_air = 50;%[50];%50;%20;%0;%50;%50;%200;% [100]

%% 150 180 280 45
% T_air_Hot = 45;%T_ambient;%270;%25;%120;%100;%25;%250;%200;%40;%

%%

Tback = 60;%T_Blanket0+5;%T_Blanket0+15;%45;%T_Blanket0;

hback = 50;%50000;%

T_ink0 = 30;%85;%75;% C

Mw_water = 18.01528;

Mw_PG = 76.09;

L_Water = 2257e3;

L_PG = 914e3;

%% Print Setup

DropVolume = 6.25;%8.0;%11.0;%7.63;%7.04;%6.3;% % pL [6.3]

SolidsPercentages = 17;%25.8;%10.1;%27.8;%28.0;%99.9999;%99.9;%28;%10;%30;%99.9;%28;%10;%30;% % % [17]

PGPercentages = 16;%8;%35;%4.5;%15.6;%0;%0;%4.5;%35;%0;%15;%0;%4.5;%35;%0;%0;%15;% % % [16]

WaterPercentages = 100 - SolidsPercentages - PGPercentages; % %

ResolutionX = 600;%600;% dpi [600]

ResolutionY = 1200; % dpi [1200]

% dt = 1e-3; % sec

%% Blanket Setup

dlmax = 10e-3; %m

k0 = 0.2; %W/m-K
ro0 = 1100; %kg/m3
Cp0 = 1500; %J/kg-K

k_Solids = k0;
ro_Solids = ro0;
Cp_Solids = Cp0;

k_Water = 0.609;
ro_Water = 1000;
Cp_Water = 4184;

k_PG = 0.147;
ro_PG = 1040;
Cp_PG = 2500;

nz = length(dz); %ul

T_Blanket = T_Blanket0*ones(nz,1); %C

k = k0*ones(nz,1); %W/m-K
ro = ro0*ones(nz,1); %kg/m3
Cp = Cp0*ones(nz,1); %J/kg-K

Ch = ro.*Cp.*dz; %J/m2-K (Capacitance per unit area)
ih = dz./k; %(m2-K)/W


% W_Ink_0 = (DropVolume*1e-15)/((0.0254/ResolutionX)*(0.0254/ResolutionY))*(AreaCoverage/100);

time = Length / Print_Speed; % s


% dz_Ink_min = (DropVolume*1e-15)/((0.0254/ResolutionX)*(0.0254/ResolutionY))*(AreaCoverage/100)*(SolidsPercentages/100);
% k_Ink_ofMin = k_PG;
% ih_Ink_min = dz_Ink_min/k_Ink_ofMin;

tauu = Ch(2:nz)  .*(ih(2:nz)/2+ih(1:nz-1)/2); %s
taud = Ch(1:nz-1).*(ih(2:nz)/2+ih(1:nz-1)/2); %s

Nyquist = 2;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tau1min = Ch(1)*(ih(1)/2+ih_Ink_min/2); %s
% taunzmin = Ch(nz)*(ih(nz)/2+1/hback); %s
% ChInkmin = ro.*Cp.*dz_Ink_min;
% tauInkmin1 = ChInkmin*(ih(1)/2+ih_Ink_min/2); %s
% tauInkmin2 = ChInkmin*(1/h_air+ih_Ink_min/2); %s
% dtmin = min([tauInkmin1 ; tauInkmin2 ; tauu ; taud ; tau1min ; taunzmin]); %s
% dt = dtmin/(Nyquist*(1+10*eps)); %s
% dl = dt * Print_Speed; %m
% Nindlmax = ceil(dlmax/dl);
% dl = dlmax / Nindlmax;
% dt = dl / Print_Speed;
% au = dt./tauu; %ul
% ad = dt./taud; %ul
% anz = dt/taunzmin; %ul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N_time = int32(time / dt);
%% Analysis

% W_Solids = ones(N_time+1,1)*(DropVolume*1e-15)/((0.0254/ResolutionX)*(0.0254/ResolutionY))*(AreaCoverage/100)*(SolidsPercentages/100);
% W_PG = zeros(N_time+1,1);
% W_Water = zeros(N_time+1,1);
% T_Ink = zeros(N_time+1,1);
% E_L = zeros(N_time+1,1);

% W_Solids(1) = (DropVolume*1e-15)/((0.0254/ResolutionX)*(0.0254/ResolutionY))*(AreaCoverage/100)*(SolidsPercentages/100);
% W_PG(1) = (DropVolume*1e-15)/((0.0254/ResolutionX)*(0.0254/ResolutionY))*(AreaCoverage/100)*(PGPercentages/100); % m
% W_Water(1) = (DropVolume*1e-15)/((0.0254/ResolutionX)*(0.0254/ResolutionY))*(AreaCoverage/100)*(WaterPercentages/100); % m
W_Solids_100 = (DropVolume*1e-15)/((0.0254/ResolutionX)*(0.0254/ResolutionY))*(SolidsPercentages/100);
W_PG_100 = (DropVolume*1e-15)/((0.0254/ResolutionX)*(0.0254/ResolutionY))*(PGPercentages/100); % m
W_Water_100 = (DropVolume*1e-15)/((0.0254/ResolutionX)*(0.0254/ResolutionY))*(WaterPercentages/100); % m

W_Solids(1) = HelpLayer;%W_Solids_100*(AreaCoverage_0/100);
W_PG(1) = 0;%W_PG_100*(AreaCoverage_0/100); % m
W_Water(1) = 0;%W_Water_100*(AreaCoverage_0/100); % m
T_Ink(1) = T_Blanket0;%T_ink0;%


for idx=1:7
W_PB(idx).Solids  = W_Solids_100*(AreaCoverage_PB(idx)/100);
W_PB(idx).PG = W_PG_100*(AreaCoverage_PB(idx)/100);
W_PB(idx).Water = W_Water_100*(AreaCoverage_PB(idx)/100);
end



% T_Boil_initial = T_Boil((PGPercentages/Mw_PG)/(PGPercentages/Mw_PG+WaterPercentages/Mw_water+eps))

X_air = Humidity_Ratio_Water(T_Room_ambient,RH_ambient); % kg/kg

Q_From_Dryer = 0;

han = waitbar(0,'0.00%','Name','Drying...  :-)');

timePass(1)=0;
% dt=dt0;
% while i<=N_time+1
i=2;
for idx = 1:7
    PB_Ready(idx) = true;
end
DemiLayer=true;

z_Wet_Initial = W_PG(1)+W_Water(1);

for idx=1:7
    z_Wet_Initial = z_Wet_Initial + W_PB(idx).PG + W_PB(idx).Water;
end

firstTimeInDryer = true;
%IRDryerON = false;

%T_Blanket_AirDryer_vector=[];
%q_Blanket_IRDryer = 0;



DoMagic(timePass, i, time, PB_Start ,W_Solids, W_PG, W_Water, k_Solids, k_PG, k_Water, Ch, ih, nz, hback, ro, Cp, h_air, tauu, taud, Nyquist, Print_Speed, dlmax, han, T_Ink, PGPercentages, WaterPercentages, v_air, X_air, Mw_PG, Mw_water, T_Blanket, Tback, Radiation_start, AirIPU_Start, AirIPU_end, h_IPU, T_IPU, ro_Solids, ro_PG, ro_Water, Cp_Solids, Cp_PG, Cp_Water, L_Water, L_PG, Q_From_Dryer, Width, AirDryer_Start, T_Press_Ambient, PB_Ready, AreaCoverage_PB, DemiLayer, W_PB, Radiation_end, firstTimeInDryer, LayerNumber, IndexOfFirstAbsorptionLayer, q_radiation_Layer, dzAbsorbing, h_air_IR, T_IR, AirDryer_end, Length, dz, FigureNumber, z_Wet_Initial, TotalCoverage, q_radiation_Blanket, AbsorbedEnergyDevidedByInputEnergy, IR_EffectiveLength, Length_IR_Dryer, LengthDryer, p0, TestNumber, T_ink0)

end

function DoMagic(timePass, i, time, PB_Start ,W_Solids, W_PG, W_Water, k_Solids, k_PG, k_Water, Ch, ih, nz, hback, ro, Cp, h_air, tauu, taud, Nyquist, Print_Speed, dlmax, han, T_Ink, PGPercentages, WaterPercentages, v_air, X_air, Mw_PG, Mw_water, T_Blanket, Tback, Radiation_start, AirIPU_Start, AirIPU_end, h_IPU, T_IPU, ro_Solids, ro_PG, ro_Water, Cp_Solids, Cp_PG, Cp_Water, L_Water, L_PG, Q_From_Dryer, Width, AirDryer_Start, T_Press_Ambient, PB_Ready, AreaCoverage_PB, DemiLayer, W_PB, Radiation_end, firstTimeInDryer, LayerNumber, IndexOfFirstAbsorptionLayer, q_radiation_Layer, dzAbsorbing, h_air_IR, T_IR, AirDryer_end, Length, dz, FigureNumber, z_Wet_Initial, TotalCoverage, q_radiation_Blanket, AbsorbedEnergyDevidedByInputEnergy, IR_EffectiveLength, Length_IR_Dryer, LengthDryer, p0, TestNumber, T_ink0)
T_Blanket_IRDryer_vector=[];
q_Blanket_AirDryer = 0;
T_InkPB(1) = T_ink0;
T_InkPB(2) = T_ink0;
T_InkPB(3) = T_ink0;
T_InkPB(4) = T_ink0;
T_InkPB(5) = T_ink0;
T_InkPB(6) = T_ink0;
T_InkPB(7) = T_ink0;

while timePass(i-1)<=time
    
    if timePass(i-1)>=PB_Start(1)*time && PB_Ready(1) && AreaCoverage_PB(1)>0
        [dz_Ink, T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(1, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, dz_Ink, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water);
    elseif timePass(i-1)>=PB_Start(2)*time && PB_Ready(2) && AreaCoverage_PB(2)>0
        [dz_Ink, T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(2, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, dz_Ink, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water);
    elseif timePass(i-1)>=PB_Start(3)*time && PB_Ready(3) && AreaCoverage_PB(3)>0
        [dz_Ink, T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(3, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, dz_Ink, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water);
    elseif timePass(i-1)>=PB_Start(4)*time && PB_Ready(4) && AreaCoverage_PB(4)>0
        [dz_Ink, T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(4, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, dz_Ink, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water);
    elseif timePass(i-1)>=PB_Start(5)*time && PB_Ready(5) && AreaCoverage_PB(5)>0
        [dz_Ink, T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(5, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, dz_Ink, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water);
    elseif timePass(i-1)>=PB_Start(6)*time && PB_Ready(6) && AreaCoverage_PB(6)>0
        [dz_Ink, T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(6, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, dz_Ink, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water);
    elseif timePass(i-1)>=PB_Start(7)*time && PB_Ready(7) && AreaCoverage_PB(7)>0
        [dz_Ink, T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(7, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, dz_Ink, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water);
    end
        
    
    dz_Ink = (W_Solids(i-1)+W_PG(i-1)+W_Water(i-1));
    k_Ink = (W_Solids(i-1)*k_Solids+W_PG(i-1)*k_PG+W_Water(i-1)*k_Water)/dz_Ink;
    ih_Ink = dz_Ink/k_Ink;
    
    tau1min = Ch(1)*(ih(1)/2+ih_Ink/2); %s
    taunzmin = Ch(nz)*(ih(nz)/2+1/hback); %s
    ChInkmin = ro.*Cp.*dz_Ink;
    tauInkmin1 = ChInkmin*(ih(1)/2+ih_Ink/2); %s
    tauInkmin2 = ChInkmin*(1/h_air+ih_Ink/2); %s
    dtmin = min([tauInkmin1 ; tauInkmin2 ; tauu ; taud ; tau1min ; taunzmin]); %s
    
    dt = dtmin/(Nyquist*(1+10*eps)); %s
    dl = dt * Print_Speed; %m
    Nindlmax = ceil(dlmax/dl);
    dl = dlmax / Nindlmax;
    dt = dl / Print_Speed;
    
    au = dt./tauu; %ul
    ad = dt./taud; %ul
    anz = dt/taunzmin; %ul
%     
    timePass(i)=timePass(i-1)+dt;
    
% for i=2:N_time+1
    
%     ratio = double(i)/double(N_time);
    ratio = timePass(i-1)/time;
    
    if mod(ratio*1000,1)/10<0.1
        
        waitbar(ratio,han,[num2str(ratio*100,'%1.1f') '%']);
        
    end
%     s3=size(T_Ink(i-1))
    Xs_PG = Humidity_Ratio_PG(T_Ink(i-1),100,PGPercentages/(PGPercentages+WaterPercentages+eps));
    
    kgp = ( 25 + 19 * (v_air+Print_Speed)) / 3600; % kg/m^2-s
    
    mdot = kgp * ( Xs_PG - X_air ); % kg/m^2-s
    
    T_Boil_out = T_Boil((PGPercentages/Mw_PG)/(PGPercentages/Mw_PG+WaterPercentages/Mw_water+eps));
    
    VaporRatio_out = VaporRatio(T_Boil_out);
    
    mdot_PG = mdot*VaporRatio_out; % kg/m^2-s
    
    mdot_water = mdot*(1-VaporRatio_out); % kg/m^2-s
    
    if W_Water(i-1)>0
        
        W_Water(i) = W_Water(i-1) - mdot_water / ro_Water * dt; % m
        
        if W_Water(i) < 0
            
            W_Water(i) = 0;
            
        end
        
    else
        
        W_Water(i) = 0;
        
    end
    
    if W_PG(i-1)>0
        
        W_PG(i) = W_PG(i-1) - mdot_PG / ro_PG * dt; % m
        
        if W_PG(i) < 0
            
            W_PG(i) = 0;
            
        end
        
    else
        
        W_PG(i) = 0;
        
    end
    
    W_Solids(i) = W_Solids(i-1);
    
    
    PGPercentages = W_PG(i)/(W_PG(i)+W_Water(i)+W_Solids(i));
    
    WaterPercentages = W_Water(i)/(W_PG(i)+W_Water(i)+W_Solids(i));
    
    
    tau1 = Ch(1).*(ih(1)/2+ih_Ink/2); %s
    a1 = dt./tau1; %ul
    dT_Blanket = diff(T_Blanket(:,i-1));
    T_Blanket(:,i) = T_Blanket(:,i-1) + [a1.*(T_Ink(i-1)-T_Blanket(1,i-1)) ; -au.*dT_Blanket] + [ad.*dT_Blanket ; anz*(Tback-T_Blanket(nz,i-1))];
    %     T_Blanket(LayerNumber,i) = T_Blanket(LayerNumber,i) + Q_radiation/(Length*Width)*dt/(ro(LayerNumber).*Cp(LayerNumber).*dz(LayerNumber));
    ih_up = 1/h_air;
%     i>=round(Radiation_start*N_time) && i<round(Radiation_end*N_time)
    if timePass(i)>=(Radiation_start*time) && timePass(i)<(Radiation_end*time)
        if firstTimeInDryer
            DryerStartIndex = i;
            firstTimeInDryer = false;
        end
%         Rad_Length = Length_IR_Dryer*Rad_Ratio/Rad_Number
        T_Blanket(LayerNumber,i) = T_Blanket(LayerNumber,i) + (q_radiation_Layer(LayerNumber-IndexOfFirstAbsorptionLayer+1)*dt/(ro(1).*Cp(1).*dzAbsorbing*1e-6))';
%         T_Ink(i) = T_Ink(i) + (Q_radiation_Ink)/(Length_IR_Dryer*Width)*dt/(ro_Ink.*Cp_Ink.*dz_Ink);
        ih_up = 1/h_air_IR;
        T_air_Hot = T_IR;
        T_Blanket_IRDryer_vector = [T_Blanket_IRDryer_vector T_Blanket(1,i)];
        T_Blanket_IRDryer_mean = mean(T_Blanket_IRDryer_vector);
        q_Blanket_IRDryer = 1/ih_up*(T_air_Hot-T_Blanket_IRDryer_mean);
        IRDryerON = true;
    elseif timePass(i)>=(AirIPU_Start(1)*time) && timePass(i)<(AirIPU_end(1)*time)
        ih_up = 1/h_IPU;
        T_air_Hot = T_IPU(1);
    elseif timePass(i)>=(AirIPU_Start(2)*time) && timePass(i)<(AirIPU_end(2)*time)
        ih_up = 1/h_IPU;
        T_air_Hot = T_IPU(2);
    elseif timePass(i)>=(AirIPU_Start(3)*time) && timePass(i)<(AirIPU_end(3)*time)
        ih_up = 1/h_IPU;
        T_air_Hot = T_IPU(3);
    elseif timePass(i)>=(AirIPU_Start(4)*time) && timePass(i)<(AirIPU_end(4)*time)
        ih_up = 1/h_IPU;
        T_air_Hot = T_IPU(4);
    elseif timePass(i)>=(AirIPU_Start(5)*time) && timePass(i)<(AirIPU_end(5)*time)
        ih_up = 1/h_IPU;
        T_air_Hot = T_IPU(5);
    elseif timePass(i)>=(AirIPU_Start(6)*time) && timePass(i)<(AirIPU_end(6)*time)
        ih_up = 1/h_IPU;
        T_air_Hot = T_IPU(6);
    elseif timePass(i)>=(AirIPU_Start(7)*time) && timePass(i)<(AirIPU_end(7)*time)
        ih_up = 1/h_IPU;
        T_air_Hot = T_IPU(7);
    elseif timePass(i)>=(AirDryer_Start*time) && timePass(i)<(AirDryer_end*time)
        if firstTimeInDryer
            DryerStartIndex = i;
            firstTimeInDryer = false;
        end
        ih_up = 1/h_AirDryer;
        T_air_Hot = T_AirDryer;
        T_Blanket_AirDryer_vector = [T_Blanket_AirDryer_vector T_Blanket(1,i)];
        T_Blanket_AirDryer_mean = mean(T_Blanket_AirDryer_vector);
        q_Blanket_AirDryer = 1/ih_up*(T_air_Hot-T_Blanket_AirDryer_mean);
    else
        T_air_Hot = T_Press_Ambient;
    end
    
    ro_Ink = (W_Solids(i-1)*ro_Solids+W_PG(i-1)*ro_PG+W_Water(i-1)*ro_Water)/dz_Ink;
    Cp_Ink = (W_Solids(i-1)*Cp_Solids+W_PG(i-1)*Cp_PG+W_Water(i-1)*Cp_Water)/dz_Ink;
    Ch_Ink = ro_Ink.*Cp_Ink.*dz_Ink;
    tau_Ink = Ch_Ink.*(ih(1)/2+ih_Ink/2);
    a_Ink = dt./tau_Ink; %ul
%     ih_up = 1/h_air;
    %         if i>=round(Air_Start*N_time) && i<round(Air_end*N_time)
    %     ih_up = 1/h_IPU;
    %     T_air_Hot = T_IPU;
    %         else
    %     T_air_Hot = T_Press_Ambient;
    %     end;
    
    %     ih_up
    %     ih_Ink
    tau_Ink_up = Ch_Ink.*(ih_up+ih_Ink/2); %s
    a_Ink_up = dt./tau_Ink_up; %ul
    
    T_Ink(i) = T_Ink(i-1) + a_Ink.*(T_Blanket(1,i-1)-T_Ink(i-1))...
        - (W_Water(i-1)-W_Water(i))*L_Water/Cp_Ink/dz_Ink...
        - (W_PG(i-1)-W_PG(i))*L_PG/Cp_Ink/dz_Ink + a_Ink_up*(T_air_Hot-T_Ink(i-1));
    
    q_L (i) = (W_Water(i-1)-W_Water(i))*ro_Water*L_Water...
        + (W_PG(i-1)-W_PG(i))*ro_PG*L_PG;
    
    Q_From_Dryer = Q_From_Dryer + h_air * (T_air_Hot-T_Ink(i-1)) * Width * Print_Speed * dt;
i=i+1;
end

delete(han);

if DemiLayer
W_Solids = W_Solids - W_Solids(1);
end

% Q_L = ((W_Water(1)-W_Water(end))*ro_Water*L_Water + (W_PG(1)-W_PG(end))*ro_PG*L_PG)
% sum(q_L)
%% Print results

% Power = Width*Print_Speed/time*(dt)^2*sum(q_L)/1000

Time_Vector_Temp = linspace(0,time,length(W_Solids));

Length_vector = Print_Speed * Time_Vector_Temp;

% figure(1);plot(Length_vector,(W_Solids+W_PG+W_Water)*1e6,'.-');
% figure(11);plot(linspace(0,time,length(W_Solids)),(W_Solids+W_PG+W_Water)*1e6,'.-');
% 
% xlim([Length_vector(1) Length_vector(end)]);
% xlim([Time_Vector_Temp(1) Time_Vector_Temp(end)]);
% ylim([0,max((W_Solids+W_PG+W_Water)*1e6*1.1)]);
% 
% xlabel('Time (sec)');
% 
% ylabel('Ink Thickness (um)');

% title('Drying');

Power = cumsum(q_L)*Width*Length/1000;
% Q = Power(end)
% figure(3);plot(linspace(0,time,length(W_Solids)),Power,'.-');
% 
% % xlim([Length_vector(1) Length_vector(end)]);
% xlim([Time_Vector_Temp(1) Time_Vector_Temp(end)]);
% ylim([0,max(Power)*1.1]);
% ylabel('Cumulative Power (kW)');
% xlabel('Time (sec)');


% figure(2);plot(Length_vector,(W_Solids+W_PG+W_Water)*1e6,'.-',Length_vector,(W_Solids+W_PG)*1e6,'.-',Length_vector,W_Solids*1e6,'.-');
%
% ylim([0,max((W_Solids+W_PG+W_Water)*1e6*1.1)]);
%
% figure(3);semilogy(Length_vector,(W_Water)*1e6,'.-');
%
% ylim([0,max((W_Solids+W_PG+W_Water)*1e6*1.1)]);

% Total_Time = double(N_time)*dt;
Depth = [0 ; cumsum(dz)];% ; sum(dz)+100e-6];%;%-dz0/2;
% Depth = Depth(end:-1:1);
% TimeVector = 0:dt:Total_Time;
TimeVector = timePass;
LengthVector = TimeVector * Print_Speed;

% Tmax = max(max(T_Blanket))
T_plot = [T_Blanket ; T_Ink];
% T_plot = T_plot;%(end:-1:1,:);
Depth_plot = -Depth*1000;%(end:-1:1)*1000;

figure(FigureNumber);
subplot(3,1,3);
% imagesc(LengthVector,Depth*1000,T_Blanket);
surf(LengthVector,Depth_plot,T_plot);
caxis([min(min(T_Blanket)) max(max(T_Blanket))]);
xlim([LengthVector(1) LengthVector(end)]);
ylim([Depth_plot(end) Depth_plot(1)]);
shading flat;
view(2);
% imagesc(T_Blanket);
xlabel('Dryer Length (m)');
ylabel('Blanket Thickness (um)');


subplot(3,1,2);
% plot(LengthVector,T_Blanket(1,:),'.-');
plot(LengthVector,T_Ink,'.-');
xlim([LengthVector(1) LengthVector(end)]);
xlabel('Dryer Length (m)');
ylabel('Ink Temperature (C)');

TinkEND = T_Ink(end);

z_Water_Initial = W_Water(1);
z_PG_Initial = W_PG(1);

for idx =1:7
    z_Water_Initial = z_Water_Initial + W_PB(idx).Water;
    z_PG_Initial = z_PG_Initial + W_PB(idx).PG;
end


WaterInLiquid = W_Water./(W_PG+W_Water)*100;
z_Wet_finals = W_PG+W_Water;

DriedLiquid = (1-z_Wet_finals/z_Wet_Initial)*TotalCoverage;
DriedWater = (1-W_Water(end)/z_Water_Initial)*100;%%
DriedPG = (1-W_PG(end)/z_PG_Initial)*100;%%
%% Power calculation
% Q_AirDryerBlanket = h_AirDryer*(T_Dryer-T_Blanket) *Width*LengthAirDryer
Q_IRDryer = IRDryerON*(q_radiation_Blanket*AbsorbedEnergyDevidedByInputEnergy(end)*IR_EffectiveLength+q_Blanket_IRDryer)*Width*Length_IR_Dryer*1e-3;
Q_AirDryer = q_Blanket_AirDryer*Width*LengthDryer*1e-3;
%%
Q_IR_Wall = p0/1000*Width*Length_IR_Dryer;
PowerDisplay=['Q_IRDryer=' num2str(Q_IRDryer,'%1.0f') 'kW ; Q_AirDryer=' num2str(Q_AirDryer,'%1.0f') 'kW ; Q_IR_Wall=' num2str(Q_IR_Wall,'%1.0f') 'kW'];
ToDisplay = ['Test Number=' num2str(TestNumber,'%1.0f') ' ; DriedWater=' num2str(DriedWater,'%1.0f') '% ; DriedPG=' num2str(DriedPG,'%1.0f') '% ; T final=' num2str(TinkEND,'%1.0f') 'C ; IPU=' num2str(DriedLiquid(DryerStartIndex),'%1.0f') '% ; Total=' num2str(DriedLiquid(end),'%1.0f') '%'];

subplot(3,1,1);plot(LengthVector,(W_Solids+W_PG+W_Water)*1e6,'.-');
% figure(31);plot(LengthVector,(W_Solids+W_PG+W_Water)*1e6,'.-');
xlim([LengthVector(1) LengthVector(end)]);
xlabel('Dryer Length (m)');
ylabel('Ink Thickness (um)');
% ylim([0 max((W_Solids+W_PG+W_Water)*1e6)]);
ylim([0 max((W_Solids+W_PG+W_Water)*1e6+eps)]);
title(ToDisplay);
%%
disp([ToDisplay ' ; ' PowerDisplay]);
%%

% figure(41);plot(LengthVector,DriedLiquid,'.-',LengthVector,WaterInLiquid,'.-');
% xlim([LengthVector(1) LengthVector(end)]);
% ylim([0 TotalCoverage]);
% xlabel('Dryer Length (m)');
% ylabel('Dried Liq. [Bl] & Water/Liq. [Gr] (%)');
% title(['Dried liquid = ' num2str(DriedLiquid(end),'%1.1f') '% ; Water in liquid = ' num2str(WaterInLiquid(1),'%1.1f') '%']);


%% Outputs
% disp(['T_final=' num2str(TinkEND,'%1.1f') 'C']);
% disp(['Test_Number=' num2str(TestNumber,'%1.0f') ' ; T_final=' num2str(TinkEND,'%1.1f') 'C ; IPU=' num2str(DriedLiquid(DryerStartIndex),'%1.0f') '% ; Total=' num2str(DriedLiquid(end),'%1.0f') '%']);
%%
% Tbup = T_Blanket(1,end);
% Tbdo = T_Blanket(end,end);
% Q_From_Dryer_kW = Q_From_Dryer/1000;
% Q_To_Blanket_kW = ro0*Width*z*Print_Speed*Cp0*(T_Ink(end)-T_Ink(1))/1000;
% Q_To_Ink_kW = ro_Ink*Width*dz_Ink*Print_Speed*Cp_Ink*(mean(T_Blanket(:,end))-mean(T_Blanket(:,1)))/1000;

% close Figure(468);

toc;

end


function [dz_Ink, T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(idx, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, dz_Ink, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water)
    if ~DemiLayer
        dz_Ink = (W_Solids(i-1)+W_PG(i-1)+W_Water(i-1));
        T_Ink(i-1) = (T_Ink(i-1)*(ro_Solids*Cp_Solids*W_Solids(i-1)+ro_PG*Cp_PG*W_PG(i-1)+ro_Water*Cp_Water*W_Water(i-1))+...
                      T_InkPB(idx)*(ro_Solids*Cp_Solids*W_PB(idx).Solids+ro_PG*Cp_PG*W_PB(idx).PG+ro_Water*Cp_Water*W_PB(idx).Water))/...
                     (ro_Solids*Cp_Solids*W_Solids(i-1)+ro_PG*Cp_PG*W_PG(i-1)+ro_Water*Cp_Water*W_Water(i-1)+...
                      ro_Solids*Cp_Solids*W_PB(idx).Solids+ro_PG*Cp_PG*W_PB(idx).PG+ro_Water*Cp_Water*W_PB(idx).Water);
        W_Solids(i-1) = W_Solids(i-1) + W_PB(idx).Solids;
        W_PG(i-1) = W_PG(i-1) + W_PB(idx).PG;
        W_Water(i-1) = W_Water(i-1) + W_PB(idx).Water;
    else
        T_Ink(i-1)=T_InkPB(idx);
        W_Solids(1:i-1) = 0;%W_Solids(i-1) + W_SolidsPB1;
        W_Solids(i-1) = W_Solids(i-1) + W_PB(idx).Solids;
        W_PG(i-1) = W_PG(i-1) + W_PB(idx).PG;
        W_Water(i-1) = W_Water(i-1) + W_PB(idx).Water;
        DemiLayer=false;
    end
    PB_Ready(idx) = false;

end

function X = Humidity_Ratio_Water(Temperature,RH)

Pa = 101325;% Pa

Pw = ((0.0506)+exp(20.386-5132/(273.15+Temperature))*133.322368); % Pa

X = 0.62198 * (RH/100)*Pw / (Pa - (RH/100)*Pw); % kg/kg

end

function X = Humidity_Ratio_PG(Tin,RH,PGinW)

Pa = 740;% Torr

T = ([32 40:20:400]-32)*(5/9);

PG = [0 40 60 70 75 80 85 90 95 96 97 98 99 100.1]/100;

i_PG = find((PGinW-PG)<0,1) - 1;

P(1,:) = [4.5 6.2 14 26 49 90 160 250 390 600 850 1 1 1 1 1 1 1 1 1];
P(2,:) = [3.9 5.5 12 23 42 78 140 220 340 520 760 1 1 1 1 1 1 1 1 1];
P(3,:) = [3.5 5 10 20 37 70 120 200 300 480 690 1000 1 1 1 1 1 1 1 1];
P(4,:) = [2.9 4.1 8.8 17 32 60 100 170 270 410 600 900 1 1 1 1 1 1 1 1];
P(5,:) = [2.6 3.6 7.7 16 29 54 90 160 240 380 550 800 1 1 1 1 1 1 1 1];
P(6,:) = [2.2 3.1 6.7 14 26 48 80 145 220 340 500 750 1 1 1 1 1 1 1 1];
P(7,:) = [1.7 2.4 5.2 11 20 38 65 120 170 270 400 600 850 1 1 1 1 1 1 1];
P(8,:) = [1 1.8 3.9 8 16 28 50 83 140 220 320 480 680 900 1 1 1 1 1 1];
P(9,:) = [1 1 2 4.3 8.4 16 28 49 78 140 190 300 440 580 810 1 1 1 1 1];
P(10,:) = [1 1 1.6 3.3 6.7 14 23 40 63 110 170 250 370 500 705 1000 1 1 1 1];
P(11,:) = [1 1 1 2.3 4.6 9 17 29 49 80 140 200 300 410 595 830 1 1 1 1];
P(12,:) = [1 1 1 1 3 6.1 12 21 37 60 97 165 230 330 490 700 960 1 1 1];
P(13,:) = [1 1 1 1 1.8 3.9 7.1 14 23 40 64 110 165 240 350 500 690 950 1 1];
P(14,:) = [1 1 1 1 1 1 2.5 5 9.1 18 29 50 82 140 200 300 420 630 900 1];

P1 = P(i_PG,:);
P2 = P(i_PG+1,:);

log10P1 = log10(P1);
log10P2 = log10(P2);

% Tin
i = find((Tin-T)<0,1) - 1;
% Tin

% s1=size((T(i+1)-Tin))
% s2=size(log10P1(i))
% s3=size(Tin-T(i))
% s4=size(log10P1(i+1))

log10P1out = ((T(i+1)-Tin)*log10P1(i)+(Tin-T(i))*log10P1(i+1))...
    /(T(i+1)-T(i));

log10P2out = ((T(i+1)-Tin)*log10P2(i)+(Tin-T(i))*log10P2(i+1))...
    /(T(i+1)-T(i));

log10Pout = ((PG(i_PG+1)-PGinW)*log10P1out+(PGinW-PG(i_PG))*log10P2out)...
    /(PG(i_PG+1)-PG(i_PG));

P1test = P1(i);
P2test = P2(i);
Pw1 = 10^log10P1out;
Pw2 = 10^log10P2out;

Pw = 10^log10Pout;

Mw_water = 18.01528;
Mw_PG = 76.09;
Mw_dryAir = 29;
% PGinW=0;
X = (Mw_water*(1-PGinW)+Mw_PG*PGinW)/Mw_dryAir * (RH/100)*Pw / (Pa - (RH/100)*Pw); % kg/kg

end

function T_Boil = T_Boil(HumectantsPercentage)

R_Liquid = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

T_Boil_Liquid = [100 103 106 108 112 117 122 131 142 158 188];

if HumectantsPercentage < R_Liquid(end)
    
    Ind_min_R_Liquid = find((HumectantsPercentage-R_Liquid)<0,1) - 1;
    
    T_Boil = ((R_Liquid(Ind_min_R_Liquid+1)-HumectantsPercentage)*T_Boil_Liquid(Ind_min_R_Liquid)...
        +(HumectantsPercentage-R_Liquid(Ind_min_R_Liquid))*T_Boil_Liquid(Ind_min_R_Liquid+1))...
        /(R_Liquid(Ind_min_R_Liquid+1)-R_Liquid(Ind_min_R_Liquid));
    
elseif HumectantsPercentage >= R_Liquid(end)
    
    T_Boil = T_Boil_Liquid(end);
    
end

end

function VaporRatio = VaporRatio(T_Boil)

T_Boil_Liquid = [100 131 145 154 161 167 171 176 180 184 188];

VaporRatio_Vapor = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

if T_Boil < T_Boil_Liquid(end)
    
    Ind_min_T_Boil_Liquid = find((T_Boil-T_Boil_Liquid)<0,1) - 1;
    
    VaporRatio = ((T_Boil_Liquid(Ind_min_T_Boil_Liquid+1)-T_Boil)*VaporRatio_Vapor(Ind_min_T_Boil_Liquid)...
        +(T_Boil-T_Boil_Liquid(Ind_min_T_Boil_Liquid))*VaporRatio_Vapor(Ind_min_T_Boil_Liquid+1))...
        /(T_Boil_Liquid(Ind_min_T_Boil_Liquid+1)-T_Boil_Liquid(Ind_min_T_Boil_Liquid));
    
elseif T_Boil == T_Boil_Liquid(end)
    
    VaporRatio = VaporRatio_Vapor(end);
    
end

end

