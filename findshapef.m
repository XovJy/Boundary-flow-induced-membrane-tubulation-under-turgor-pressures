clear;

fp = 0;          % pressure
ka = 1;     
ri = 1e-3;
rb = 1.0;         % base radius
mu = 1;           % visocisty
Ga1 = 0;          % normal friction
Ga2 = 0.01;       % tangential friction
A0 = 0.01;


% velocity
vu = -1.5;

% simulation parameters
RTol = 1e-6;
Nmax = 5000;
dt = 0.1;
T = 100;

% save path
filename = ['./SolData/dynamic_vu_',num2str(vu),...
            '_Fp_',num2str(fp),...
	        '_T_',num2str(T),...
	        '_dt_',num2str(dt),...
            '_date_',date,'.mat'];

data = load("./InitData/init_vu0_0_tarvu_-10_Ga2_0.01.mat");
index = -vu * 10 + 1;
sol = data.sol_all{index, 2};

sol_all = cell(3000,2);
Tr = linspace(0,T,3000);
j = 0;

for t = 0:dt:T
    t
    disp(sol.y(5,end));
    disp(sol.y(9,end));
    uc = sol.x;
    r0c = sol.y(3,:);
    z0c = sol.y(4,:);
    psi0c = sol.y(1,:);
    h0 = sol.parameters;
    area = sol.y(5,end);

    yeq = @(u,y,para)shapef(u,y,para,fp,ka,mu,Ga1,Ga2,uc,r0c,z0c,psi0c,h0,dt);
    ybc = @(ya,yb,para)twobcf(ya,yb,para,vu,ka,rb,ri);
    jach = @(u,y,para) jac(u,y,para,fp,ka,mu,Ga1,Ga2,uc,r0c,z0c,psi0c,h0,dt);
    opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax,'FJacobian',jach,'Vectorized','on');
    sol = bvp5c(yeq,ybc,sol,opts);
   
    if min(abs(Tr-t))<0.5*dt
        j = j+1;
        sol_all{j,1} = t;
        sol_all{j,2} = sol;
    end
    save(filename,'sol_all');
end


