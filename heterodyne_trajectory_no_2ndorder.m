function [m2x, m2y, m2z]=heterodyne_trajectory_no_2ndorder(Jxeff,Jyeff,Jzeff,gamma, timearray, dt, initial_state, Nlength, factormeasured)
%-------------------------------------------------------------------------%
%   Time evolve a quantum trajectory with heterodyne unravelling using the
%   truncated cumulant hierarchy trajectory method up to second order. This
%   time evolution does not include the Wiener noise terms in the second 
%   order cumulant equations (see arXiv:2209.13377v3 for more information).
%   
%Parameters:
%   Jxeff           Coupling strength in the x-direction
%   Jyeff           Coupling strength in the y-direction
%   Jzeff           Coupling strength in the z-direction
%   gamma           Spontaneous local dissipation strength
%   timearray       Time values at which to save results
%   dt              Numerical time step
%   initial_state   Initial state of the simulation
%   Nlength         Lattice dimension
%   factormeasured  Measurement efficiency
%
%Returns:
%   m2x             Structure factor (k=0) in the x-direction
%   m2y             Structure factor (k=0) in the y-direction
%   m2z             Structure factor (k=0) in the z-direction
%-------------------------------------------------------------------------%

sqdt=sqrt(dt);                                                              % Square root of numerical timestep: 
gammam = factormeasured*gamma;                                              % gamma_measured (fraction that is measured properly)
Nspins=Nlength^2;                                                           % Total number of spins in the lattice
NBmat=sparse(Nspins,Nspins);                                                % Connection matrix of coupled/ interacting spins

% Set the coupled/neighbouring spins for each spin in a 2D square lattice:
for n=1:Nspins
    leftn=n-1; if mod(leftn,Nlength)==0,leftn=leftn+Nlength; end
    rightn=n+1;if mod(rightn,Nlength)==1, rightn=rightn-Nlength;end
    downn=n+Nlength; if (mod(downn,Nspins)<=Nlength && mod(downn,Nspins)~=0),downn=downn-Nspins;end 
    upn=n-Nlength; if (mod(upn,Nspins)>Nspins-Nlength||mod(upn,Nspins)==0),upn=upn+Nspins; end    
    NBmat(n,leftn)=1;NBmat(n,rightn)=1;NBmat(n,downn)=1;NBmat(n,upn)=1;
end
%NOTE: in this example we work with a 2D square checkerboard lattice. One
%can however opt for any lattice geometry and dimension.

%set initial state
sigx=initial_state(:,1)+0.01*(randn(Nspins,1));                             %Add noise (as a precaution) to avoid hitting an unstable solution.
sigy=initial_state(:,2)+0.01*(randn(Nspins,1));                             %Add noise (as a precaution) to avoid hitting an unstable solution.
sigz=zeros(Nspins,1)+0.01*(randn(Nspins,1));                                %Add noise (as a precaution) to avoid hitting an unstable solution.

% Normalize the magnetizations to make sure the spin states are 
% on the Bloch sphere:
normalise = sqrt(sigx.^2 + sigy.^2 + sigz.^2);
sigx = sigx./normalise;
sigy = sigy./normalise;
sigz = sigz./normalise;

% Set the values for the cumulants:
fxfx=zeros(Nspins,Nspins);
fyfy=zeros(Nspins,Nspins);
fzfz=zeros(Nspins,Nspins);
fxfy=zeros(Nspins,Nspins);
fyfz=zeros(Nspins,Nspins); 
fzfx=zeros(Nspins,Nspins); 

% Set the storage arrays for the magnetization and store 
% the intitial state at t = 0
sigxres=NaN(Nspins,length(timearray));sigxres(:,1)=sigx;
sigyres=NaN(Nspins,length(timearray));sigyres(:,1)=sigy;
sigzres=NaN(Nspins,length(timearray)); sigzres(:,1)=sigz; savequantumcorrs=true;
 if(savequantumcorrs) %put to false to restrict the memory that the output requires
fxfxresultarray=NaN(Nspins,Nspins,length(timearray)); fxfxresultarray(:,:,1)=fxfx;
fyfyresultarray=NaN(Nspins,Nspins,length(timearray)); fyfyresultarray(:,:,1)=fyfy;
fzfzresultarray=NaN(Nspins,Nspins,length(timearray)); fzfzresultarray(:,:,1)=fzfz;
fxfyresultarray=NaN(Nspins,Nspins,length(timearray)); fxfyresultarray(:,:,1)=fxfy;
fyfzresultarray=NaN(Nspins,Nspins,length(timearray)); fyfzresultarray(:,:,1)=fyfz;
fzfxresultarray=NaN(Nspins,Nspins,length(timearray)); fzfxresultarray(:,:,1)=fzfx;
 end

% Start the time evolution
time=timearray(1);
while time<timearray(end)
    
    % From this point onwards the various contributions to the equations of
    % motion are calculated. We refer to arXiv:2209.13377v3 for the 
    % equations. 
    sxsy=fxfy+sigx.*sigy';                                                  % Calculate the (x,y)-correlations
    sysz=fyfz+sigy.*sigz';                                                  % Calculate the (y,z)-correlations
    szsx=fzfx+sigz.*sigx';                                                  % Calculate the (z,x)-correlations                                                  % Calculate the (x,y)-correlations


    dWx=sqdt*randn(Nspins,1); dWy=sqdt*randn(Nspins,1);                     % The random Wiener noise

    % Deterministic and stochastic contributions to the equations for the 
    % magnetization in the x,y and z-direction.
    d_sigx_det=-0.5*gamma*sigx+2*Jyeff*diag(NBmat*sysz)-2*Jzeff*diag(sysz*NBmat); 
    d_sigx_stoch=sqrt(gammam/2)*(fxfx+diag(1+sigz-sigx.^2))*dWx-sqrt(gammam/2)*(fxfy-diag(sigx.*sigy))*dWy;
    d_sigy_det=-0.5*gamma*sigy+2*Jzeff*diag(NBmat*szsx)-2*Jxeff*diag(szsx*NBmat);
    d_sigy_stoch=sqrt(gammam/2)*(fxfy'-diag(sigx.*sigy))*dWx-sqrt(gammam/2)*(diag(1+sigz-sigy.^2)+fyfy)*dWy;
    d_sigz_det=-gamma*(sigz+1)+2*Jxeff*diag(NBmat*sxsy)-2*Jyeff*diag(sxsy*NBmat);
    d_sigz_stoch=sqrt(gammam/2)*(fzfx-diag(sigx.*(1+sigz)))*dWx+sqrt(gammam/2)*(diag(sigy.*(1+sigz))-fyfz')*dWy;
    
    % Contributions to the equations for the second order cumulants:

    %evolution disregarding s'=m, m'=s situations  
    d_fxfx_nonnb=-gamma*fxfx+2*Jyeff*((NBmat*(sigy)).*fzfx+sigz.*(NBmat*(fxfy'))+fzfx'.*((sigy')*NBmat)+((fxfy)*NBmat).*sigz')-2*Jzeff*((NBmat*(sigz)).*fxfy'+sigy.*(NBmat*(fzfx))+fxfy.*(sigz'*NBmat)+(fzfx'*NBmat).*sigy');
    d_fyfy_nonnb=-gamma*fyfy+2*Jzeff*((NBmat*(sigz)).*fxfy+sigx.*(NBmat*(fyfz'))+fxfy'.*((sigz')*NBmat)+((fyfz)*NBmat).*sigx')-2*Jxeff*((NBmat*(sigx)).*fyfz'+sigz.*(NBmat*(fxfy))+fyfz.*(sigx'*NBmat)+(fxfy'*NBmat).*sigz');
    d_fzfz_nonnb=-2*gamma*(fzfz)+2*Jxeff*((NBmat*(sigx)).*fyfz+sigy.*(NBmat*(fzfx'))+fyfz'.*((sigx')*NBmat)+((fzfx)*NBmat).*sigy')-2*Jyeff*((NBmat*(sigy)).*fzfx'+sigx.*(NBmat*(fyfz))+fzfx.*(sigy'*NBmat)+(fyfz'*NBmat).*sigx');
    d_fxfy_nonnb=-gamma*fxfy-2*Jxeff*((sigx'*NBmat).*fzfx'+(fxfx*NBmat).*sigz')+2*Jyeff*((NBmat*(sigy)).*fyfz'+sigz.*(NBmat*(fyfy)))+2*Jzeff*(-(NBmat*(sigz)).*fyfy-sigy.*(NBmat*(fyfz'))+fxfx.*(sigz'*NBmat)+(fzfx'*NBmat).*sigx');
    d_fyfz_nonnb=-1.5*gamma*fyfz-2*Jyeff*((sigy'*NBmat).*fxfy'+(fyfy*NBmat).*sigx')+2*Jzeff*((NBmat*(sigz)).*fzfx'+sigx.*(NBmat*(fzfz)))+2*Jxeff*(-(NBmat*(sigx)).*fzfz-sigz.*(NBmat*(fzfx'))+fyfy.*(sigx'*NBmat)+(fxfy'*NBmat).*sigy');
    d_fzfx_nonnb=-1.5*gamma*fzfx-2*Jzeff*((sigz'*NBmat).*fyfz'+(fzfz*NBmat).*sigy')+2*Jxeff*((NBmat*(sigx)).*fxfy'+sigy.*(NBmat*(fxfx)))+2*Jyeff*(-(NBmat*(sigy)).*fxfx-sigx.*(NBmat*(fxfy'))+fzfz.*(sigy'*NBmat)+(fyfz'*NBmat).*sigz');

    %correcting for s'=m,m'=s situations.  
    d_fxfx_corr=2*Jyeff*NBmat.*(-fyfz'.*sigx'-fzfx.*sigy'-sigx.*fyfz-sigy.*fzfx'-sigz.*sigy'.*sigx'-sigz'.*sigy.*sigx)-2*Jzeff*NBmat.*(-fxfy'.*sigz'-fyfz.*sigx'-sigx.*fyfz'-sigz.*fxfy-sigy.*sigz'.*sigx'-sigx.*sigz.*sigy'); %Correction to account for s'=m, s=m' cases
    d_fyfy_corr=2*Jzeff*NBmat.*(-fzfx'.*sigy'-fxfy.*sigz'-sigy.*fzfx-sigz.*fxfy'-sigx.*sigz'.*sigy'-sigx'.*sigz.*sigy)-2*Jxeff*NBmat.*(-fyfz'.*sigx'-fzfx.*sigy'-sigy.*fzfx'-sigx.*fyfz-sigz.*sigx'.*sigy'-sigy.*sigx.*sigz');
    d_fzfz_corr=2*Jxeff*NBmat.*(-fxfy'.*sigz'-fyfz.*sigx'-sigz.*fxfy-sigx.*fyfz'-sigy.*sigx'.*sigz'-sigy'.*sigx.*sigz)-2*Jyeff*NBmat.*(-fzfx'.*sigy'-fxfy.*sigz'-sigz.*fxfy'-sigy.*fzfx-sigx.*sigy'.*sigz'-sigz.*sigy.*sigx');
    d_fxfy_corr=-2*Jxeff*NBmat.*(-2*sigx.*fzfx'+(1-sigx.^2).*sigz')+2*Jyeff*NBmat.*(-2*fyfz'.*sigy'+sigz.*(1-sigy.^2)')+2*Jzeff*NBmat.*(fyfz.*sigy'+fyfy.*sigz'-sigz.*fxfx-sigx.*fzfx+sigy.*sigy'.*sigz'-sigx.*sigx'.*sigz);
    d_fyfz_corr=-2*Jyeff*NBmat.*(-2*sigy.*fxfy'+(1-sigy.^2).*sigx')+2*Jzeff*NBmat.*(-2*fzfx'.*sigz'+sigx.*(1-sigz.^2)')+2*Jxeff*NBmat.*(fzfx.*sigz'+fzfz.*sigx'-sigx.*fyfy-sigy.*fxfy+sigz.*sigz'.*sigx'-sigy.*sigy'.*sigx);
    d_fzfx_corr=-2*Jzeff*NBmat.*(-2*sigz.*fyfz'+(1-sigz.^2).*sigy')+2*Jxeff*NBmat.*(-2*fxfy'.*sigx'+sigy.*(1-sigx.^2)')+2*Jyeff*NBmat.*(fxfy.*sigx'+fxfx.*sigy'-sigy.*fzfz-sigz.*fyfz+sigx.*sigx'.*sigy'-sigz.*sigz'.*sigy);
    
    % Ito contributions to the equations of the second order cumulants:
    d_fxfx_ito=-0.5*gamma*(fxfx+diag(1+sigz-sigx.^2))*(fxfx+diag(1+sigz-sigx.^2))'-0.5*gamma*(fxfy-diag(sigx.*sigy))*(fxfy-diag(sigx.*sigy))';
    d_fyfy_ito=-0.5*gamma*(fxfy'-diag(sigx.*sigy))*(fxfy'-diag(sigx.*sigy))'-0.5*gamma*(diag(1+sigz-sigy.^2)+fyfy)*(diag(1+sigz-sigy.^2)+fyfy)';
    d_fzfz_ito=-0.5*gamma*(fzfx-diag(sigx.*(1+sigz)))*(fzfx-diag(sigx.*(1+sigz)))'-0.5*gamma*(diag(sigy.*(1+sigz))-fyfz')*(diag(sigy.*(1+sigz))-fyfz')';
    d_fxfy_ito=-0.5*gamma*(fxfx+diag(1+sigz-sigx.^2))*(fxfy'-diag(sigx.*sigy))'-0.5*gamma*(fxfy-diag(sigx.*sigy))*(diag(1+sigz-sigy.^2)+fyfy)';
    d_fyfz_ito=-0.5*gamma*(fxfy'-diag(sigx.*sigy))*(fzfx-diag(sigx.*(1+sigz)))'+0.5*gamma*(diag(1+sigz-sigy.^2)+fyfy)*(diag(sigy.*(1+sigz))-fyfz')';
    d_fzfx_ito=-0.5*gamma*(fzfx-diag(sigx.*(1+sigz)))*(fxfx+diag(1+sigz-sigx.^2))'+0.5*gamma*(diag(sigy.*(1+sigz))-fyfz')*(fxfy-diag(sigx.*sigy))';
        
    % Update the values of the magnetization and the cumulants:
    sigx=sigx+d_sigx_det*dt+d_sigx_stoch;
    sigy=sigy+d_sigy_det*dt+d_sigy_stoch;
    sigz=sigz+d_sigz_det*dt+d_sigz_stoch;
    fxfx=fxfx+d_fxfx_nonnb*dt+d_fxfx_corr*dt+d_fxfx_ito*dt; fxfx=fxfx-diag(diag(fxfx));
    fyfy=fyfy+d_fyfy_nonnb*dt+d_fyfy_corr*dt+d_fyfy_ito*dt; fyfy=fyfy-diag(diag(fyfy));
    fzfz=fzfz+d_fzfz_nonnb*dt+d_fzfz_corr*dt+d_fzfz_ito*dt; fzfz=fzfz-diag(diag(fzfz));
    fxfy=fxfy+d_fxfy_nonnb*dt+d_fxfy_corr*dt+d_fxfy_ito*dt; fxfy=fxfy-diag(diag(fxfy));
    fyfz=fyfz+d_fyfz_nonnb*dt+d_fyfz_corr*dt+d_fyfz_ito*dt; fyfz=fyfz-diag(diag(fyfz));
    fzfx=fzfx+d_fzfx_nonnb*dt+d_fzfx_corr*dt+d_fzfx_ito*dt; fzfx=fzfx-diag(diag(fzfx));
    
    % Increment time:
    time=time+dt;

    %-----------------------------------
    %NOTE: You can opt to include checks on the physicallity of the 
    %expectation values. In case the values become unphysical,
    %you can correct them by constraining the value. This tends 
    %to allow the trajectory to remain stable.

    %An example is a check of the norm:
    nor = sqrt(sigx.*sigx + sigy.*sigy + sigz.*sigz);
    
    %renormalise sigx, sigy and sigz if their norm is larger than a
    %threshold, e.g. n_t = 1.1. 
    n_t = 1.1;
    i_t = nor > n_t;
    
    sigx(i_t) = sigx(i_t)./nor(i_t);
    sigy(i_t) = sigy(i_t)./nor(i_t);
    sigz(i_t) = sigz(i_t)./nor(i_t);
    %-----------------------------------
    

    % Save the results of the simulation at the predetermined times:
    step=find(abs(timearray-time)<dt/2); 
    if ~isempty(step)
        try
           assert(all(isfinite(sigz)),'Trajectory unstable. Evolution aborted');
        catch
            time = timearray(end);
            
        end
        sigxres(:,step)=sigx;
        sigyres(:,step)=sigy;
        sigzres(:,step)=sigz;
        
        if(savequantumcorrs)
        fxfxresultarray(:,:,step)=fxfx;
        fyfyresultarray(:,:,step)=fyfy;
        fzfzresultarray(:,:,step)=fzfz;
        fxfyresultarray(:,:,step)=fxfy;
        fyfzresultarray(:,:,step)=fyfz;
        fzfxresultarray(:,:,step)=fzfx;
        end
                    
    end
    
end

% NOTE: In this example the results will be preprocessed. Usually, this 
% will be necessary for larger lattices as the memory used to store/save 
% the data can become very large very quickly. 
% In case of very large lattices one can opt to preprocess the results at 
% each timestep and not save them for every individual spin and spin-spin 
% correlation. You can change the functionality to what is optimal for your
% calculations.

% PREPROCESSING:

m2x = zeros(1, length(timearray));
m2y = zeros(1, length(timearray));
m2z = zeros(1, length(timearray));

% Calculate observables averaged over the lattice sites for each time step:
for tt=1:length(timearray)

    % Sum of classical correlations (disregarding on-site correlations):
    Xcorrmatclassical=sigxres(:,tt).*(sigxres(:,tt).');                     % In the x-direction
    m2xclassical=sum(sum(Xcorrmatclassical-diag(diag(Xcorrmatclassical)))); % In the x-direction

    % Sum of quantum correlations (disregarding on-site correlations):
    Xcorrmatquantum=squeeze(fxfxresultarray(:,:,tt));                       % In the x-direction
    m2xquantum=sum(sum(Xcorrmatquantum-diag(diag(Xcorrmatquantum))));       % In the x-direction

    % Structure factor (k=0) in the x-direction:
    m2x(tt)=(m2xclassical+m2xquantum + Nspins)/Nspins^2;                    


    % Sum of classical correlations (disregarding on-site correlations):
    Ycorrmatclassical=sigyres(:,tt).*(sigyres(:,tt).');                     % In the y-direction
    m2yclassical=sum(sum(Ycorrmatclassical-diag(diag(Ycorrmatclassical)))); % In the y-direction

    % Sum of quantum correlations (disregarding on-site correlations):
    Ycorrmatquantum=squeeze(fyfyresultarray(:,:,tt));                       % In the y-direction
    m2yquantum=sum(sum(Ycorrmatquantum-diag(diag(Ycorrmatquantum))));       % In the y-direction

    % Structure factor (k=0) in the y-direction:
    m2y(tt)=(m2yclassical+m2yquantum + Nspins)/Nspins^2;

    % Sum of classical correlations (disregarding on-site correlations):
    Zcorrmatclassical=sigzres(:,tt).*(sigzres(:,tt).');                     % In the z-direction
    m2zclassical=sum(sum(Zcorrmatclassical-diag(diag(Zcorrmatclassical)))); % In the z-direction

    % Sum of quantum correlations (disregarding on-site correlations):
    Zcorrmatquantum=squeeze(fzfzresultarray(:,:,tt));                       % In the z-direction
    m2zquantum=sum(sum(Zcorrmatquantum-diag(diag(Zcorrmatquantum))));       % In the z-direction

    % Structure factor (k=0) in the z-direction:
    m2z(tt)=(m2zclassical+m2zquantum + Nspins)/Nspins^2;

end

%NOTE: From the above it is already clear that one can also extract
%correlations based on their nature (i.e. classical or quantum).
%Furthermore, one can also extract other observables from the remaining
%cumulants that are calculated via the hierarchy of equations, yielding
%structure factors (and others) not oriented only along a direction x, 
% y or z, but e.g. along any orientation in the (x,y)-plain.

end
