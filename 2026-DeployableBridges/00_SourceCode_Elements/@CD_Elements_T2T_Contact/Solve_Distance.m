function distance=Solve_Distance(obj,xSet1,xSet2)

    tri1=reshape(xSet1',[1,9]);
    tri2=reshape(xSet2',[1,9]);

    [distance,P1,P2]=obj.simdTriTri2(tri1,tri2);

end