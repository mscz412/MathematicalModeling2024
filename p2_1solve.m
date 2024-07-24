LB=0;
UB=100000;
options=gaoptimset('populationtype','doublevector');
options=gaoptimset(options,'populationsize',100);
options=gaoptimset(options,'PlotFcns',@gaplotbestf);
options=gaoptimset(options,'generations',200);
options=gaoptimset(options,'stallgenlimit',inf);
[x,fval]=ga(@p2_1,1,[],[],[],[],LB,UB,[],options);
