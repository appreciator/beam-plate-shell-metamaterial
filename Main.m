function Main()
% An ODE-driven Level-Set Density Method for Topology Optimization of Shell Structures
% Zuyu Li, Email:lizuyu0091@163.com; Yang Liu, Email:285548029@qq.com
%% Parameter Setting
modelPara = struct('Model',1,'penal',3,'nLoop',200,'nVolRelax',50,'dt',0.1,'reg',0.1,'nReg',4);
%% Preprocess
[modelFEM_init,modelLS_init] = Model(modelPara);
modelResponse_init = FEMSolver(modelFEM_init,modelLS_init,modelPara);
%% Optimization
LSDSolver(modelFEM_init,modelLS_init,modelResponse_init,modelPara);
end