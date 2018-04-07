function B = trainTree(Pos,Nex)
x = [Pos;Nex];
y = [ones(size(Pos,1),1); -1*ones(size(Nex,1),1)];
B = TreeBagger(10,x,y,'Method','regression');

end