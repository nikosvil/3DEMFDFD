function rel_er_10_nodes_ba(dfs,hx,hy,hz,Ex_xcoor,Ex_ycoor,Ey_xcoor,Ey_ycoor,zcoor,p_rETx,p_rETy,p_rEx_Total,p_rEy_Total)

fprintf('Coordinates | Ex_Total | ETx | Real relError\n')
fprintf('10 nodes on the left\n')
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx-3*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs-1),p_rETx(dfs-1),100*abs(p_rETx(dfs-1)-p_rEx_Total(dfs-1))/abs(p_rETx(dfs-1)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx-5*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs-2),p_rETx(dfs-2),100*abs(p_rETx(dfs-2)-p_rEx_Total(dfs-2))/abs(p_rETx(dfs-2)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx-7*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs-3),p_rETx(dfs-3),100*abs(p_rETx(dfs-3)-p_rEx_Total(dfs-3))/abs(p_rETx(dfs-3)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx-9*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs-4),p_rETx(dfs-4),100*abs(p_rETx(dfs-4)-p_rEx_Total(dfs-4))/abs(p_rETx(dfs-4)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx-11*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs-5),p_rETx(dfs-5),100*abs(p_rETx(dfs-5)-p_rEx_Total(dfs-5))/abs(p_rETx(dfs-5)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx-13*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs-6),p_rETx(dfs-6),100*abs(p_rETx(dfs-6)-p_rEx_Total(dfs-6))/abs(p_rETx(dfs-6)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx-15*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs-7),p_rETx(dfs-7),100*abs(p_rETx(dfs-7)-p_rEx_Total(dfs-7))/abs(p_rETx(dfs-7)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx-17*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs-8),p_rETx(dfs-8),100*abs(p_rETx(dfs-8)-p_rEx_Total(dfs-8))/abs(p_rETx(dfs-8)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx-19*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs-9),p_rETx(dfs-9),100*abs(p_rETx(dfs-9)-p_rEx_Total(dfs-9))/abs(p_rETx(dfs-9)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx-21*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs-10),p_rETx(dfs-10),100*abs(p_rETx(dfs-10)-p_rEx_Total(dfs-10))/abs(p_rETx(dfs-10)))
fprintf('10 nodes on the right\n')
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx+3*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs+2),p_rETx(dfs+2),100*abs(p_rETx(dfs+2)-p_rEx_Total(dfs+2))/abs(p_rETx(dfs+2)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx+5*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs+3),p_rETx(dfs+3),100*abs(p_rETx(dfs+3)-p_rEx_Total(dfs+3))/abs(p_rETx(dfs+3)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx+7*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs+4),p_rETx(dfs+4),100*abs(p_rETx(dfs+4)-p_rEx_Total(dfs+4))/abs(p_rETx(dfs+4)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx+9*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs+5),p_rETx(dfs+5),100*abs(p_rETx(dfs+5)-p_rEx_Total(dfs+5))/abs(p_rETx(dfs+5)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx+11*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs+6),p_rETx(dfs+6),100*abs(p_rETx(dfs+6)-p_rEx_Total(dfs+6))/abs(p_rETx(dfs+6)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx+13*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs+7),p_rETx(dfs+7),100*abs(p_rETx(dfs+7)-p_rEx_Total(dfs+7))/abs(p_rETx(dfs+7)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx+15*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs+8),p_rETx(dfs+8),100*abs(p_rETx(dfs+8)-p_rEx_Total(dfs+8))/abs(p_rETx(dfs+8)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx+17*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs+9),p_rETx(dfs+9),100*abs(p_rETx(dfs+9)-p_rEx_Total(dfs+9))/abs(p_rETx(dfs+9)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx+19*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs+10),p_rETx(dfs+10),100*abs(p_rETx(dfs+10)-p_rEx_Total(dfs+10))/abs(p_rETx(dfs+10)))
fprintf('%4.1f %d %d %e %e %4.2f\n',Ex_xcoor*hx+21*hx/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,...
    p_rEx_Total(dfs+11),p_rETx(dfs+11),100*abs(p_rETx(dfs+11)-p_rEx_Total(dfs+11))/abs(p_rETx(dfs+11)))

fprintf('****************************************\n')
fprintf('****************************************\n')
fprintf('Coordinates | Ey_Total | ETy | relError\n')
fprintf('10 nodes on the left\n')
% EY nodes to the left
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx-2*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs-2),p_rETy(dfs-2),100*abs(p_rETy(dfs-2)-p_rEy_Total(dfs-2))/abs(p_rETy(dfs-2)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx-3*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs-3),p_rETy(dfs-3),100*abs(p_rETy(dfs-3)-p_rEy_Total(dfs-3))/abs(p_rETy(dfs-3)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx-4*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs-4),p_rETy(dfs-4),100*abs(p_rETy(dfs-4)-p_rEy_Total(dfs-4))/abs(p_rETy(dfs-4)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx-5*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs-5),p_rETy(dfs-5),100*abs(p_rETy(dfs-5)-p_rEy_Total(dfs-5))/abs(p_rETy(dfs-5)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx-6*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs-6),p_rETy(dfs-6),100*abs(p_rETy(dfs-6)-p_rEy_Total(dfs-6))/abs(p_rETy(dfs-6)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx-7*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs-7),p_rETy(dfs-7),100*abs(p_rETy(dfs-7)-p_rEy_Total(dfs-7))/abs(p_rETy(dfs-7)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx-8*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs-8),p_rETy(dfs-8),100*abs(p_rETy(dfs-8)-p_rEy_Total(dfs-8))/abs(p_rETy(dfs-8)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx-9*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs-9),p_rETy(dfs-9),100*abs(p_rETy(dfs-9)-p_rEy_Total(dfs-9))/abs(p_rETy(dfs-9)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx-10*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs-10),p_rETy(dfs-10),100*abs(p_rETy(dfs-10)-p_rEy_Total(dfs-10))/abs(p_rETy(dfs-10)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx-11*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs-11),p_rETy(dfs-11),100*abs(p_rETy(dfs-11)-p_rEy_Total(dfs-11))/abs(p_rETy(dfs-11)))
fprintf('10 nodes on the right\n')
% EY nodes to the right
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx+2*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs+2),p_rETy(dfs+2),100*abs(p_rETy(dfs+2)-p_rEy_Total(dfs+2))/abs(p_rETy(dfs+2)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx+3*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs+3),p_rETy(dfs+3),100*abs(p_rETy(dfs+3)-p_rEy_Total(dfs+3))/abs(p_rETy(dfs+3)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx+4*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs+4),p_rETy(dfs+4),100*abs(p_rETy(dfs+4)-p_rEy_Total(dfs+4))/abs(p_rETy(dfs+4)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx+5*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs+5),p_rETy(dfs+5),100*abs(p_rETy(dfs+5)-p_rEy_Total(dfs+5))/abs(p_rETy(dfs+5)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx+6*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs+6),p_rETy(dfs+6),100*abs(p_rETy(dfs+6)-p_rEy_Total(dfs+6))/abs(p_rETy(dfs+6)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx+7*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs+7),p_rETy(dfs+7),100*abs(p_rETy(dfs+7)-p_rEy_Total(dfs+7))/abs(p_rETy(dfs+7)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx+8*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs+8),p_rETy(dfs+8),100*abs(p_rETy(dfs+8)-p_rEy_Total(dfs+8))/abs(p_rETy(dfs+8)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx+9*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs+9),p_rETy(dfs+9),100*abs(p_rETy(dfs+9)-p_rEy_Total(dfs+9))/abs(p_rETy(dfs+9)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx+10*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs+10),p_rETy(dfs+10),100*abs(p_rETy(dfs+10)-p_rEy_Total(dfs+10))/abs(p_rETy(dfs+10)))
fprintf('%d %d %d %e %e %4.2f\n',Ey_xcoor*hx+11*hx, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,...
    p_rEy_Total(dfs+11),p_rETy(dfs+11),100*abs(p_rETy(dfs+11)-p_rEy_Total(dfs+11))/abs(p_rETy(dfs+11)))

end