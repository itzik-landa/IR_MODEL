function DryingSimIR2022

%ink is transparent here

%%Blanket structure

PercentFilling = 62;%%[7 20 500] %here affect wave

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


AreaCoverage_PB(1) = 0;
AreaCoverage_PB(2) = 100;
AreaCoverage_PB(3) = 0;
AreaCoverage_PB(4) = 0;%100;%[0 100]
AreaCoverage_PB(5) = 50;%[50 50]
AreaCoverage_PB(6) = 0;%50;%[100 100]
AreaCoverage_PB(7) = 0;


HelpLayer = 10e-6; %if I run without ink, allows simulation

NumberOfIPUs = 7;

T_IPU_All = 100;%70;%90;%165;%133;%200;%    %IPU like in hot air, usualy 100
T_Blanket0 = 68;%70;%80;%                   %initial condition, currently open loop. 68 about reality
Print_Speed_SPH = 6500;%13000;%10000;%
T_IR = 50;%25;% [120];                      %IR air temerature. not to change.                      
h_air_IR = 50;                              %convection coefficient between air and blanket. not to change.    
h_IPU = 200;                                %convection coefficient between air and IPU. not to change.    

%% Other Dryer Setup
%% 0 1.5
v_air = 0;%-Print_Speed;%1.5;%0;%1.5; % m/s   %air speed 
RH_ambient = 50; % %
T_Room_ambient = 25; % %
T_Press_Ambient = 45;%25;% [45];
h_air = 50;%[50];%50;%20;%0;%50;%50;%200;% [100]
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

%% Blanket Params
dlmax = 10e-3; %m

k0 = 0.2; %W/m-K
ro0 = 1100; %kg/m3
Cp0 = 1500; %J/kg-K

k_Water = 0.609;
ro_Water = 1000;
Cp_Water = 4184;

k_PG = 0.147;
ro_PG = 1040;
Cp_PG = 2500;

%%
LengthDryer = 2.0;

PostDryer = 0.02;


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

tic;

%% Dryer Length

TotalCoverage = sum(AreaCoverage_PB);

%% Test Matrix

dzRelease = ReleaseThickness/ReleaseLayersNumber;
dzAbsorbing = AbsorbingLayerThickness/AbsorbingLayersNumber;
dzAboveFabric = AboveFabricLayerThickness/AboveFabricLayersNumber;
dzFabricLayer = FabricLayerThickness/FabricLayersNumber;
dzGripLayers = GripLayerThickness/GripLayersNumber;

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

p0_coeff = 2.4; %2.4
Rr = 0.05; %0.05

p0 = p0_coeff*39e3/LengthDryer; %power of the lamp per length. The coefficent can be adjusted. Look how much in the new
q_radiation_Blanket = GetRadiation(p0, Rr); 

params = ['p0-coeff=' num2str(p0_coeff) ', Rr=' num2str(Rr)];

for idx = 1:7
    PB_Start(idx) = (0.375*(idx-1)+0.2+0.175/2)/Length;
    T_IPU (idx) = T_IPU_All;
    AirIPU_Start(idx) = (0.375*(idx-1)+0.2)/Length;
    AirIPU_end(idx) = (0.375*(idx-1)+0.2)/Length;
end

AirDryer_Start = 0.375*NumberOfIPUs/Length;
AirDryer_end = (0.375*NumberOfIPUs+LengthDryer)/Length*0;

Length_IR_Dryer = LengthDryer;
Width = 1;%m
Area_IR_Dryer = Length_IR_Dryer * Width;
Area_Dryer = LengthDryer * Width;

Radiation_start = 0.375*NumberOfIPUs/Length;
Radiation_end = (0.375*NumberOfIPUs+Length_IR_Dryer)/Length;
Rad_Ratio = 3/4;
IR_EffectiveLength = Rad_Ratio;
IndexOfFirstAbsorptionLayer = ReleaseLayersNumber + 1;%#
Alpha=0.276;%1/(%*u)[0.276]
for j=1:AbsorbingLayersNumber+1
    AbsorptionDepth = (j-1) * dzAbsorbing;%um
    AbsorbedEnergyDevidedByInputEnergy(j) = 1 - exp(-Alpha*PercentFilling/100*AbsorptionDepth);%%
end
AbsorptionPerLayer = diff(AbsorbedEnergyDevidedByInputEnergy);%%
LayerNumber = IndexOfFirstAbsorptionLayer:(IndexOfFirstAbsorptionLayer+AbsorbingLayersNumber-1);%
q_radiation_Layer = q_radiation_Blanket*AbsorptionPerLayer;
figure(468);plot(linspace(0,AbsorbingLayerThickness,length(AbsorbedEnergyDevidedByInputEnergy)),AbsorbedEnergyDevidedByInputEnergy,'.-');








%% Blanket Setup
k_Solids = k0;
ro_Solids = ro0;
Cp_Solids = Cp0;

nz = length(dz); %ul

T_Blanket = T_Blanket0*ones(nz,1); %C

k = k0*ones(nz,1); %W/m-K
ro = ro0*ones(nz,1); %kg/m3
Cp = Cp0*ones(nz,1); %J/kg-K

Ch = ro.*Cp.*dz; %J/m2-K (Capacitance per unit area)
ih = dz./k; %(m2-K)/W

time = Length / Print_Speed; % s


tauu = Ch(2:nz)  .*(ih(2:nz)/2+ih(1:nz-1)/2); %s
taud = Ch(1:nz-1).*(ih(2:nz)/2+ih(1:nz-1)/2); %s

Nyquist = 2;

%% Analysis

W_Solids_100 = (DropVolume*1e-15)/((0.0254/ResolutionX)*(0.0254/ResolutionY))*(SolidsPercentages/100);
W_PG_100 = (DropVolume*1e-15)/((0.0254/ResolutionX)*(0.0254/ResolutionY))*(PGPercentages/100); % m
W_Water_100 = (DropVolume*1e-15)/((0.0254/ResolutionX)*(0.0254/ResolutionY))*(WaterPercentages/100); % m

W_Solids(1) = HelpLayer;%W_Solids_100*(AreaCoverage_0/100);
W_PG(1) = 0;%W_PG_100*(AreaCoverage_0/100); % m
W_Water(1) = 0;%W_Water_100*(AreaCoverage_0/100); % m
T_Ink(1) = T_Blanket0;%T_ink0;%

z_Wet_Initial = W_PG(1)+W_Water(1);

for idx=1:7
    W_PB(idx).Solids  = W_Solids_100*(AreaCoverage_PB(idx)/100);
    W_PB(idx).PG = W_PG_100*(AreaCoverage_PB(idx)/100);
    W_PB(idx).Water = W_Water_100*(AreaCoverage_PB(idx)/100);
    PB_Ready(idx) = true;
    z_Wet_Initial = z_Wet_Initial + W_PB(idx).PG + W_PB(idx).Water;
end

X_air = Humidity_Ratio_Water(T_Room_ambient,RH_ambient); % kg/kg

han = waitbar(0,'0.00%','Name','Drying...  :-)');



DoMagic(time, PB_Start ,W_Solids, W_PG, W_Water, k_Solids, k_PG, k_Water, Ch, ih, nz, hback, ro, Cp, h_air, tauu, taud, Nyquist, Print_Speed, dlmax, han, T_Ink, PGPercentages, WaterPercentages, v_air, X_air, Mw_PG, Mw_water, T_Blanket, Tback, Radiation_start, AirIPU_Start, AirIPU_end, h_IPU, T_IPU, ro_Solids, ro_PG, ro_Water, Cp_Solids, Cp_PG, Cp_Water, L_Water, L_PG, Area_IR_Dryer, Area_Dryer, AirDryer_Start, T_Press_Ambient, PB_Ready, AreaCoverage_PB, W_PB, Radiation_end, LayerNumber, IndexOfFirstAbsorptionLayer, q_radiation_Layer, dzAbsorbing, h_air_IR, T_IR, AirDryer_end, dz, z_Wet_Initial, TotalCoverage, q_radiation_Blanket, AbsorbedEnergyDevidedByInputEnergy, IR_EffectiveLength, p0, T_ink0, params)

end

function DoMagic(time, PB_Start ,W_Solids, W_PG, W_Water, k_Solids, k_PG, k_Water, Ch, ih, nz, hback, ro, Cp, h_air, tauu, taud, Nyquist, Print_Speed, dlmax, han, T_Ink, PGPercentages, WaterPercentages, v_air, X_air, Mw_PG, Mw_water, T_Blanket, Tback, Radiation_start, AirIPU_Start, AirIPU_end, h_IPU, T_IPU, ro_Solids, ro_PG, ro_Water, Cp_Solids, Cp_PG, Cp_Water, L_Water, L_PG, Area_IR_Dryer, Area_Dryer, AirDryer_Start, T_Press_Ambient, PB_Ready, AreaCoverage_PB, W_PB, Radiation_end, LayerNumber, IndexOfFirstAbsorptionLayer, q_radiation_Layer, dzAbsorbing, h_air_IR, T_IR, AirDryer_end, dz, z_Wet_Initial, TotalCoverage, q_radiation_Blanket, AbsorbedEnergyDevidedByInputEnergy, IR_EffectiveLength, p0, T_ink0, params)
i=2;
DemiLayer=true;
firstTimeInDryer = true;
timePass(1)=0;

T_Blanket_IRDryer_vector=[];
q_Blanket_AirDryer = 0;

for idx = 1:7
    T_InkPB(idx) = T_ink0;
end

while timePass(i-1)<=time
    
    if timePass(i-1)>=PB_Start(1)*time && PB_Ready(1) && AreaCoverage_PB(1)>0
        [T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(1, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water);
    elseif timePass(i-1)>=PB_Start(2)*time && PB_Ready(2) && AreaCoverage_PB(2)>0
        [T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(2, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water);
    elseif timePass(i-1)>=PB_Start(3)*time && PB_Ready(3) && AreaCoverage_PB(3)>0
        [T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(3, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water);
    elseif timePass(i-1)>=PB_Start(4)*time && PB_Ready(4) && AreaCoverage_PB(4)>0
        [T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(4, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water);
    elseif timePass(i-1)>=PB_Start(5)*time && PB_Ready(5) && AreaCoverage_PB(5)>0
        [T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(5, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water);
    elseif timePass(i-1)>=PB_Start(6)*time && PB_Ready(6) && AreaCoverage_PB(6)>0
        [T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(6, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water);
    elseif timePass(i-1)>=PB_Start(7)*time && PB_Ready(7) && AreaCoverage_PB(7)>0
        [T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(7, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water);
    end
        
    dz_Ink = (W_Solids(i-1)+W_PG(i-1)+W_Water(i-1));
    k_Ink = (W_Solids(i-1)*k_Solids+W_PG(i-1)*k_PG+W_Water(i-1)*k_Water)/dz_Ink;
    ih_Ink = dz_Ink/k_Ink;
    
    tau1min = Ch(1)*(ih(1)/2+ih_Ink/2); 
    taunzmin = Ch(nz)*(ih(nz)/2+1/hback); 
    ChInkmin = ro.*Cp.*dz_Ink;
    tauInkmin1 = ChInkmin*(ih(1)/2+ih_Ink/2); 
    tauInkmin2 = ChInkmin*(1/h_air+ih_Ink/2); 
    dtmin = min([tauInkmin1 ; tauInkmin2 ; tauu ; taud ; tau1min ; taunzmin]); %s
    
    dt = dtmin/(Nyquist*(1+10*eps)); 
    dl = dt * Print_Speed; 
    Nindlmax = ceil(dlmax/dl);
    dl = dlmax / Nindlmax;
    dt = dl / Print_Speed;
    
    au = dt./tauu; 
    ad = dt./taud; 
    anz = dt/taunzmin;  
    timePass(i)=timePass(i-1)+dt;

    ratio = timePass(i-1)/time;
    
    if mod(ratio*1000,1)/10<0.1
        
        waitbar(ratio,han,[num2str(ratio*100,'%1.1f') '%']);
        
    end
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
    
    
    tau1 = Ch(1).*(ih(1)/2+ih_Ink/2);
    a1 = dt./tau1;
    dT_Blanket = diff(T_Blanket(:,i-1));
    T_Blanket(:,i) = T_Blanket(:,i-1) + [a1.*(T_Ink(i-1)-T_Blanket(1,i-1)) ; -au.*dT_Blanket] + [ad.*dT_Blanket ; anz*(Tback-T_Blanket(nz,i-1))];
    ih_up = 1/h_air;
    if timePass(i)>=(Radiation_start*time) && timePass(i)<(Radiation_end*time)
        if firstTimeInDryer
            DryerStartIndex = i;
            firstTimeInDryer = false;
        end
        T_Blanket(LayerNumber,i) = T_Blanket(LayerNumber,i) + (q_radiation_Layer(LayerNumber-IndexOfFirstAbsorptionLayer+1)*dt/(ro(1).*Cp(1).*dzAbsorbing*1e-6))';
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
    a_Ink = dt./tau_Ink;

    tau_Ink_up = Ch_Ink.*(ih_up+ih_Ink/2);
    a_Ink_up = dt./tau_Ink_up;
    
    T_Ink(i) = T_Ink(i-1) + a_Ink.*(T_Blanket(1,i-1)-T_Ink(i-1))...
        - (W_Water(i-1)-W_Water(i))*L_Water/Cp_Ink/dz_Ink...
        - (W_PG(i-1)-W_PG(i))*L_PG/Cp_Ink/dz_Ink + a_Ink_up*(T_air_Hot-T_Ink(i-1));

i=i+1;
end

delete(han);

if DemiLayer
W_Solids = W_Solids - W_Solids(1);
end

PrintItOut(dz, timePass, Print_Speed, T_Blanket, T_Ink, W_Water, W_PG, W_PB, z_Wet_Initial, TotalCoverage, IRDryerON, q_radiation_Blanket, AbsorbedEnergyDevidedByInputEnergy...
    , IR_EffectiveLength, q_Blanket_IRDryer,Area_Dryer, Area_IR_Dryer, q_Blanket_AirDryer, p0, DryerStartIndex, W_Solids, params);

end

function [pb_tr] = GetRadiation(p0, Rr)

%new lamp can control p0 (power - lower) and Tr (how much transmitted to
%blanket - lower)

%look for spectrum of ink. Can calculate actual Tr. Check what happens if
%only ink absorbs. maybe Tr is 0?

%p0 כמה עוצמה הגיע לפני השטח של השמיכה
%כמה בלנקט מקבל מבחינת צפיפות אנרגיה
    Te = 0.60;

    Tr = 1 - Rr;
    Tb_tr = Te * Tr; %how much transmitted to the absorbing layer. The ink could take part of it.
    pb_tr = Tb_tr * p0;
end

function PrintItOut(dz, timePass, Print_Speed, T_Blanket, T_Ink, W_Water, W_PG, W_PB, z_Wet_Initial, TotalCoverage, IRDryerON, q_radiation_Blanket, AbsorbedEnergyDevidedByInputEnergy,...
    IR_EffectiveLength, q_Blanket_IRDryer,Area_Dryer, Area_IR_Dryer, q_Blanket_AirDryer, p0, DryerStartIndex, W_Solids, params)
fig_out = 1;
TestNumber = fig_out;
FigureNumber = TestNumber;
%% Print results

Depth = [0 ; cumsum(dz)];
TimeVector = timePass;
LengthVector = TimeVector * Print_Speed;

T_plot = [T_Blanket ; T_Ink];
Depth_plot = -Depth*1000;

figure(FigureNumber);
subplot(3,1,3);
surf(LengthVector,Depth_plot,T_plot);
caxis([min(min(T_Blanket)) max(max(T_Blanket))]);
xlim([LengthVector(1) LengthVector(end)]);
ylim([Depth_plot(end) Depth_plot(1)]);
shading flat;
view(2);
xlabel('Dryer Length (m)');
ylabel('Blanket Thickness (um)');


subplot(3,1,2);
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


z_Wet_finals = W_PG+W_Water;

DriedLiquid = (1-z_Wet_finals/z_Wet_Initial)*TotalCoverage;
DriedWater = (1-W_Water(end)/z_Water_Initial)*100;%%
DriedPG = (1-W_PG(end)/z_PG_Initial)*100;%%
%% Power calculation

Q_IRDryer = IRDryerON*(q_radiation_Blanket*AbsorbedEnergyDevidedByInputEnergy(end)*IR_EffectiveLength+q_Blanket_IRDryer)*Area_IR_Dryer*1e-3;
Q_AirDryer = q_Blanket_AirDryer*Area_Dryer*1e-3;

Q_IR_Wall = p0/1000*Area_IR_Dryer;
PowerDisplay=['Q_IRDryer=' num2str(Q_IRDryer,'%1.0f') 'kW ; Q_AirDryer=' num2str(Q_AirDryer,'%1.0f') 'kW ; Q_IR_Wall=' num2str(Q_IR_Wall,'%1.0f') 'kW'];
ToDisplay = ['Test Number=' num2str(TestNumber,'%1.0f') ' ; DriedWater=' num2str(DriedWater,'%1.0f') '% ; DriedPG=' num2str(DriedPG,'%1.0f') '% ; T final=' num2str(TinkEND,'%1.0f') 'C ; IPU=' num2str(DriedLiquid(DryerStartIndex),'%1.0f') '% ; Total=' num2str(DriedLiquid(end),'%1.0f') '% ;' newline params];

subplot(3,1,1);plot(LengthVector,(W_Solids+W_PG+W_Water)*1e6,'.-');
xlim([LengthVector(1) LengthVector(end)]);
xlabel('Dryer Length (m)');
ylabel('Ink Thickness (um)');
ylim([0 max((W_Solids+W_PG+W_Water)*1e6+eps)]);
title(ToDisplay);

disp([ToDisplay ' ; ' PowerDisplay]);


toc;
end


function [T_Ink, W_Solids, W_PG, W_Water, DemiLayer, PB_Ready] = SetValues(idx, i, DemiLayer, W_Solids, W_PG, W_Water, T_Ink, T_InkPB, W_PB, PB_Ready, ro_Solids, Cp_Solids, ro_PG, Cp_PG, ro_Water, Cp_Water)
    if ~DemiLayer
        %dz_Ink = (W_Solids(i-1)+W_PG(i-1)+W_Water(i-1));
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


log10P1out = ((T(i+1)-Tin)*log10P1(i)+(Tin-T(i))*log10P1(i+1))...
    /(T(i+1)-T(i));

log10P2out = ((T(i+1)-Tin)*log10P2(i)+(Tin-T(i))*log10P2(i+1))...
    /(T(i+1)-T(i));

log10Pout = ((PG(i_PG+1)-PGinW)*log10P1out+(PGinW-PG(i_PG))*log10P2out)...
    /(PG(i_PG+1)-PG(i_PG));


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

