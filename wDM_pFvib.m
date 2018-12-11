
clear

temp=0:25:300;
press = 0;

bz_sampled = 1;

%%%%%%
%Simple Fvib
%%%%%%
z=thermalExpansion('glycine_alpha.freq','e-el-b86bpbeXDM.dat',temp,press,bz_sampled,'results')


function out = evaluateDoubleMurnaghan(x,Vm,Gm,comp,expan)
  %Murnaghan
  %F_M = @(a,x) a(1) + (a(2) * x)/a(3) .* (((a(4)./x).^(a(3)))/(a(3)nano gl - 1) + 1 ) - ((a(2) * a(4) )/(a(3) - 1));
  F_M = @(a,x) x(:,3) + (a(1) * x(:,1))/a(2) .* (((x(:,2)./x(:,1)).^(a(2)))/(a(2) - 1) + 1 ) - ((a(1) * x(:,2) )/(a(2) - 1));
  out = zeros(length(x),1);
  %length(x)
  %x(:)
  test = zeros(length(x),3);
  test(:,1) = x(:);
  test(:,2) = Vm;
  test(:,3) = Gm;
  x = test;
  %size(x,1)
  if(size(x,1)>1)
    %fprintf('Correctly found a vector\n');
    %Evaluate each element of the array
    
    for i = 1:1:size(x,1)
        %x(i,:)
        %x(:,i)
        %x(i,1)
        %x(i,2)
        %x(i,3)
      if(x(i,1)<=Vm)
        out(i) = F_M(comp,x(i,:));
      else
        out(i) = F_M(expan,x(i,:));
      end
    end
  else
    if(x<=Vm)
      out = F_M(comp,x);
    else
      out = F_M(expan,x);
    end
  end
  %return out
end

function y = thermalExpansion(freq_file,el_energy_file,Temp,Press,count_structures,dirName)

%Script to compute Free Energy at given Temp and Volumes
%Self-generates Grueneisen parameters and Fvib calcs
%clear

feedBack = 'Beginning Job';

%%%%
% First Check that the mandatory parameters were 
% initialized and set some defaults
%%%%
if(isempty(freq_file))
  feedBack = char(feedBack,'Warning! Frequency file was not set! Exiting...');
  y=feedBack;
  return;
end

if(isempty(el_energy_file))
  feedBack = char(feedBack,'Warning! Energy file was not set! Exiting...');
  y=feedBack;
  return;
end

if(isempty(Temp))
  feedBack = char(feedBack,'Warning! Temperature was not set! Assuming 0');
  Temp = 0;
end

if(isempty(Press))
  feedBack = char(feedBack,'Warning! Pressure was not set! Assuming 0');
  Press = 0;
end

if(isempty(count_structures))
  feedBack = char(feedBack,'Warning! count_structures was not set! Assuming 1');
  count_structures = 1;
end

if(isempty(dirName))
  feedBack = char(feedBack,'Warning! Directory was not set! Assuming current directory');
  dirName = pwd;
end

currentDir=pwd;
if ~strcmp(currentDir,dirName) & ~exist(dirName,'dir')
  mkdir(dirName);
end
phononDir=sprintf('%s/%s', dirName, 'phonons');
eosDir=sprintf('%s/%s', dirName, 'eosFits');
fvibDir=sprintf('%s/%s', dirName, 'fvib');
freeEnergyDir=sprintf('%s/%s', dirName, 'freeEnergy');
summaryDir=sprintf('%s/%s', dirName, 'summary');

if ~exist(phononDir, 'dir') 
    mkdir(phononDir);
end
if ~exist(eosDir, 'dir') 
    mkdir(eosDir);
end
if ~exist(fvibDir, 'dir') 
    mkdir(fvibDir);
end
if ~exist(freeEnergyDir, 'dir') 
    mkdir(freeEnergyDir);
end
if ~exist(summaryDir, 'dir') 
    mkdir(summaryDir);
end

%%%%
% Now set some variables that will be used throughout the function
%%%%

test = fopen(freq_file);
freqs = [];
tmp_freqs = [];
freq_volumes = {};
x_old = 0;
eosSet = 0;
eosType = '';
acceptedEOS = '';
options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',1000);

%constants needed for the calculation
h= 6.62607363e-34; %Js
Na= 6.0221367e23; %1/mol
kb= 1.3806488e-23; %J/K
R= kb*Na;  %J/(mol K)
c= 2.99792458e8; %m/s
%count_structures = 0;
hartreeTokJperMol = 2625.500; %kJ/mol / Hartee

%
%Step 0: Extract Fvib data from file
%
%test = fopen(freq_file);
%fvib_data = importdata(freq_file);
%fvib_volumes = fvib_data.data(:,1)';
%fvib = fvib_data.data(:,2)' .* hartreeTokJperMol;

% First the raw E(V) data points
%scatter(fvib_volumes,fvib); 
%hold on; % keep plot active while we add curves
count_freqs = 0;

tline = fgetl(test);
while ~feof(test)
    %disp(tline);
    if contains(tline,'volume')
        blah = strsplit(tline);
        vol = str2double(blah(2));
        freq_volumes = [freq_volumes vol];
        if(count_freqs > 0)
          %convert from a cell array to a regular array
          %tmp_freqs = cell2mat(tmp_freqs);
          
          %get rid of NaN in the array
          tmp_freqs(isnan(tmp_freqs)) = [];
          
          %size_tmp1 = size(tmp_freqs)
          %size_tmp2 = size(tmp_freqs2)
          %Add freq set to the total Freq count
          freqs = [freqs tmp_freqs'];
          %tmp_freqs2 = cell2mat(tmp_freqs2);
          %freqs = cat(2,freqs,tmp_freqs);
        end
        tmp_freqs = [];
        tline = fgetl(test);
        count_freqs = count_freqs + 1;
    end
    
    x = str2double(tline);
    
    tmp_freqs = [tmp_freqs x];
    
    tline = fgetl(test);
end

          
%get rid of NaN in the array
tmp_freqs(isnan(tmp_freqs)) = [];

%Add the final freq set to the total Freq count
freqs = [freqs tmp_freqs'];
%tmp_freqs2 = cell2mat(tmp_freqs2)

fclose(test);
%count_structures = count_structures / 3; 
%count_structures = 1;
%fprintf('count_structures: %i\n',count_structures);
%tmpstr = sprintf('Number of Bruellion zones sampled: %i',count_structures);
%feedBack = char(feedBack,tmpstr);

%tmp_freqs2(1,:)
%num=size(tmp_freqs2)
%dim=size(tmp_freqs2{1,1})
%freqs = zeroes(dim,num)
%for length = 1:num
%  freqs(,length) = tmp_freqs2{dim,length};
%end

freq_volumes = cell2mat(freq_volumes);
%freqs = cell2mat(freqs)


%freqs(1,:)
%size(freqs)
%Now ensure there are no imaginary frequencies
%If so, zero those and the cooresponding modes in 
%the reference frequency lists

index=find(freqs(1:end,1:end)<0);

if(~isempty(index));
  feedBack = char(feedBack,'Warning! Found an imaginary reference frequency. Zeroing...');
  freqs(index) = 0;
end




%%%%
%Now sort the frequencies acording to the volumes
%from lowest volume to greatest
%%%%

orig_freq_volumes = freq_volumes;
orig_freqs = freqs;

freq_volumes = sort(freq_volumes);
for length = 1:size(freq_volumes,2)
    pos = find(orig_freq_volumes == freq_volumes(length));
    freqs(:,length) = orig_freqs(:,pos);
end

%Now scale the frequencies
%units: m/s * (cm/m) * (1/cm) = 1/s = Hz
freqs = freqs .* ( c * 100 );

%
%Step 2: Extract electronic energy data from file.
%

%e_el = importdata('e-el-dft-alpha.dat');
e_el = importdata(el_energy_file);
%Volumes = e_el(:,1)';
%Energies = e_el(:,2)';
Volumes = e_el.data(:,1)';
Energies = e_el.data(:,2)';
[minEnergy,pos] = min(Energies);
minVolume = Volumes(pos);

%
%Step 3: Now establish the Volumes, Pressures, and Temperatures we wish to
%calculate over
%

%Volume to extrapolate from
%v_min = Energy_fit(4)-100;
%v_max = Energy_fit(4)+100;
%v_min = 300;
%v_max = 800;
%vol_array = v_min:v_max;

%
%Step 3: Now establish the Volumes, Pressures, and Temperatures we wish to
%calculate over
%


%Volume to extrapolate from
v_min = min(Volumes(:));
v_max = max(Volumes(:));
vol_array=Volumes;
dim=size(freqs,2);

%Pressure in GPa
for press = Press
  holdName = sprintf('final_resultsP%1.2f.txt', press);
  results = fopen(holdName,'w');
  fprintf(results,'#Temperature (K)     Optimal Volume (Ang^3)    Free Energy (kJ/mol)    Entropy (J/(mol K))    Enthalpy (kJ/mol)    Fvib (kJ/mol)    Internal Energy (kJ/mol)\n');
  
  %Loop over a set of temperatures
  for temp = Temp
    %Reset Fvib values to an empty array
    fprintf('\n\nTemp = %f\n',temp);
    fvib = zeros(1,dim);
    hvib = zeros(1,dim);
    svib = zeros(1,dim);
    entropy_enthalpy = zeros(1,dim);
    s_vib = zeros(1,dim);
    h_vib_term2 = zeros(1,dim);
    %cv = zeros(1,dim);
    %pv = {};
    
    %Calculate zero point energy contribution
    zpe = (sum(freqs) * h) / 2;

    %If relevant calculate other term
    if temp ~= 0
     %units: J/K * K = J
     h_vib_term2 = (h .* freqs) ./ ( exp( ( h .* freqs )./( kb * temp )) - 1);
     h_vib_term2 (~isfinite(h_vib_term2))=0;
     h_vib_term2 = sum(h_vib_term2);
     s_vib = ((h .* freqs) ./(temp .* ( exp( ( h .* freqs )./( kb * temp )) - 1))) - (( kb ) .* log( 1 - exp( ( -1 * h .* freqs )./( kb * temp ))));
     s_vib(~isfinite(s_vib))=0;
     s_vib= sum(s_vib);   
     %Cv = (Na/(count_structures*kb)).*( ( ((h .* freqs) ./(temp .* ( exp( ( h .* freqs )./( kb * temp )) - 1))).^2).* exp( (h .* freqs )./( kb * temp )));
     %Cv(~isfinite(Cv))=0;
     %Cv = sum(Cv);     
     %fprintf('entr_enth: %f kJ/mol\n',entropy_enthalpy);       
    %else
      %entropy_enthalpy = 0.0;
      %s_vib = 0;
      %h_vib_term2 = 0;
      %Cv = 0;
    end
    
    %Sum to create Fvib
    %units: (/mol) * (J + J) * (kJ / J) = kJ/mol
    svib =  Na .* (s_vib) ./ (1000 * count_structures);
    hvib =  Na .* (zpe + h_vib_term2) ./ (1000 * count_structures);
    fvib = hvib - (temp .* svib);
    
    %Perform some fits to back out values
    fvib_polyfit = polyfit(freq_volumes,fvib,3);
    poly_fvib_func = @(v) polyval(fvib_polyfit,v);
    
    hvib_polyfit = polyfit(freq_volumes,hvib,3);
    poly_hvib_func = @(v) polyval(hvib_polyfit,v);  
  
    svib_polyfit = polyfit(freq_volumes,svib,3);
    poly_svib_func = @(v) polyval(svib_polyfit,v); 
    %poly_fvib_func(freq_volumes)
    
    % First the raw E(V) data points
    %%%%
    %scatter(freq_volumes,fvib); 
    %hold on; % keep plot active while we add curves
    %plot(freq_volumes,poly_fvib_func(freq_volumes));
    %hold off;
    %%%%
    
    pv = press * Volumes * (1.0e-24) * Na;
    gibbs = poly_fvib_func(Volumes) + Energies + pv;
    
    gibbs_polyfit = polyfit(vol_array,gibbs,4);
    gibbs_polyfit_func = @(v) polyval(gibbs_polyfit,v);
    %gibbs_spline = spline(vol_array,gibbs);
    %spline_gibbs_func = @(v) ppval(gibbs_spline,v);


    [minGibbs,pos] = min(gibbs);
    if (vol_array(pos) == vol_array(end));
      pos = pos - 1;
      minVol = vol_array(pos);
      minGibbs = gibbs(pos);
    elseif (vol_array(pos) == vol_array(1));
      pos = pos + 1;
      minVol = vol_array(pos);
      minGibbs = gibbs(pos);
    else
      minVol = vol_array(pos);
    end
    
    [Vmin, Gmin, info, output] = fminbnd(gibbs_polyfit_func, v_min, v_max);

    %Now construct the double Murnaghan EOS
    %Murnaghan
    %F_M = @(a,x) a(1) + (a(2) * x)/a(3) .* (((a(4)./x).^(a(3)))/(a(3) - 1) + 1 ) - ((a(2) * a(4) )/(a(3) - 1));
    %Murnaghan with restricted volumes?
    %a(4) = vo & a(1) = e0
    %Now x(:,1) = volumes
    %x(:,2) = v0 x(:,3) = e0
    %a'(1) = a(2) a'(2) = a(3)
    F_M = @(a,x) x(:,3) + (a(1) * x(:,1))/a(2) .* (((x(:,2)./x(:,1)).^(a(2)))/(a(2) - 1) + 1 ) - ((a(1) * x(:,2) )/(a(2) - 1));
    
    vol_comp = vol_array(pos-1:end);
    vol_therm = vol_array(1:pos+1);
    gibbs_comp = gibbs(pos-1:end);
    gibbs_therm = gibbs(1:pos+1);

    w_comp = vol_comp;
    for i = 1:1:size(w_comp,2);
      if(abs(w_comp(i) - minVol) == 0)
        w_comp(i) = 1.0;
      elseif(abs(w_comp(i) - minVol) < 50)
        w_comp(i) = 0.25;
      else
        w_comp(i);
      end
    end
    w_therm = vol_therm;
    for i = 1:1:size(w_therm,2);
      if(abs(w_therm(i) - minVol) == 0)
        w_therm(i) = 1.0;
      elseif(abs(w_therm(i) - minVol) < 50)
        w_therm(i) = 0.25;
      else
        w_therm(i) = 0.1;
      end
    end
    % Set initial guess parameters and perform least squares fit on the E(V)
    % data.
    a0 = [5.0, 8.0];
    comp = zeros(size(vol_comp,2),3);
    comp(:,1) = vol_comp;
    comp(:,2) = Vmin;
    comp(:,3) = Gmin;
    therm = zeros(size(vol_therm,2),3);
    therm(:,1) = vol_therm;
    therm(:,2) = Vmin;
    therm(:,3) = Gmin;
    %comp
    %therm
    total_compression_branch = fitnlm(comp,gibbs_comp,F_M,a0,'Weights',w_comp);
    total_expansion_branch = fitnlm(therm,gibbs_therm,F_M,a0,'Weights',w_therm);

    
    compression_branch = total_compression_branch.Coefficients.Estimate;
    expansion_branch = total_expansion_branch.Coefficients.Estimate;
    
    %Free_energy = @(v) F(Energy_fit,v);
    %Now search over all available volumes for the optimal Gibbs' value
    [Vmin, Gmin, info, output] = fminbnd(@(v) evaluateDoubleMurnaghan(v,Vmin,Gmin,compression_branch,expansion_branch), v_min, v_max);
    %[Vmin, Gmin, info, output] = fminbnd(poly_gibbs_func, v_min, v_max);
    %[Vmin, Gmin, info, output] = fminbnd(poly_gibbs_func, v_min, v_max);
    
    % First the raw E(V) data points
    %scatter(Volumes,gibbs); 
    %hold on; % keep plot active while we add curves
    %plot(Volumes,spline_gibbs_func(Volumes));
    %plot(vol_comp,F_M(compression_branch,comp));
    %plot(vol_therm,F_M(expansion_branch,therm));
    %plot(vol_array,polyval(gibbs_poly,vol_array));
    %plot(vol_array,poly_gibbs_func(vol_array));
    %plot(vol_array,evaluateDoubleMurnaghan(vol_array,Vmin,Gmin,compression_branch,expansion_branch));
    %scatter(Vmin,Gmin); 
    %xlabel('Volume (Ang^3)');
    %ylabel('Energy (kJ/mol)');
    %legend('Raw Data','compression', 'expansion');
    %legend('Raw Data','compression', 'expansion', 'Piece-wise');
    %legend('Raw Data','Piece-wise','Min');
    %hold off;
    
    % Print out the optimal volume and free energy
    fprintf('Optimal Volume and Free Energy:\n%f Ang^3\t%f kJ/mol\n',Vmin,Gmin);
    Smin = poly_svib_func(Vmin);
    Hmin = poly_hvib_func(Vmin);
    Fvibmin = poly_fvib_func(Vmin);
    %CVmin = spline_cv_func(Vmin);
    PVmin = press * Vmin * (1.0e-24) * Na;
    EelMin = Gmin - PVmin - Fvibmin;
    %fprintf('Optimal Entropy, Enthalpy and Cv:\n%f J/mol\t%f kJ/mol\t%f J/(mol K)\n',Smin*1000,Hmin,CVmin);
    fprintf('Optimal Entropy, Enthalpy:\n%f J/mol\t%f kJ/mol\t\n',Smin*1000,Hmin);
    fprintf('PV term :\t%f kJ/mol,\tE_el = \t%f\n',PVmin,EelMin);
  
    %
    %Step 5: Store results in 'logical' files
    %
  
    name = sprintf( 'freeEnergyT%iP%1.2f.txt', temp, press );
    fileID = fopen(name,'w');
    fprintf(fileID,'#Min Vol: %f Ang^3  Min Gibbs: %f (kJ/mol)\n',Vmin,Gmin);
    fprintf(fileID,'#Volume (Ang^3)   Free Energy (kJ/mol)    Error (kJ/mol)\n');
    for j = 1:size(vol_array,2)
      i = vol_array(j);
      fprintf(fileID,'%f %f %f\n',i,evaluateDoubleMurnaghan(i,Vmin,Gmin,compression_branch,expansion_branch),evaluateDoubleMurnaghan(i,Vmin,Gmin,compression_branch,expansion_branch)-gibbs(j));
      %fprintf('%f %f\n',i,Free_energy(i));
    end

    fclose(fileID);
    movefile(name,freeEnergyDir);
    
    name = sprintf( 'fvibEnergyT%iP%1.2f.txt', temp, press );
    fileID = fopen(name,'w');
    fprintf(fileID,'#Volume (Ang^3)   Helmholtz Free Energy (kJ/mol)\n');
    for i = vol_array
      fprintf(fileID,'%f %f\n',i,poly_fvib_func(i));
    end
    fclose(fileID);
    movefile(name,fvibDir);
    
    %fprintf(results,"%f   %f   %f   %f   %f   %f     %f\n",temp,Vmin,Fmin,Smin*1000,Hmin,CVmin,F(Energy_fit,Vmin));
    fprintf(results,"%f   %f   %f   %f   %f   %f     %f\n",temp,Vmin,Gmin,Smin*1000,Hmin,Fvibmin,EelMin);
    
    %wMin =  ( c * 100 ) .* ref_freqs .* (( Vmin / ref_vol ).^( -1 .* grun ));

    %namePhonon = sprintf( 'phononsT%iP%1.2fV%1.2f.txt', temp, press, Vmin );
    %phonon = fopen(namePhonon,'w');
    %fprintf(phonon,'Phonons at Temperature %i K Pressure %1.2f GPa and Volume  %1.2f Ang^3\n', temp, press, Vmin );
    %for i = 1:1:size(wMin(:))
    %  fprintf(phonon,'%f\n',wMin(i));
    %end
    %fclose(phonon);
    %movefile(namePhonon,phononDir);

  end
  
  fclose(results);
  movefile(holdName,summaryDir);
end

disp('Made it to the end')
feedBack = char(feedBack,'Job Complete!');
y = feedBack;
return
end

