function [p_dev_er_rEx,p_dev_er_rEy,p_dev_er_iEx,p_dev_er_iEy,p_rel_er_rEx,p_rel_er_iEx,p_rel_er_rEy,p_rel_er_iEy]=...
    comp_dev_err(p_rEx_Total,p_rETx,p_iEx_Total,p_iETx,p_rEy_Total,p_rETy,p_iEy_Total,p_iETy)

p_dev_er_rEx=abs(p_rEx_Total-p_rETx);
p_dev_er_rEy=abs(p_rEy_Total-p_rETy);
p_dev_er_iEx=abs(p_iEx_Total-p_iETx);
p_dev_er_iEy=abs(p_iEy_Total-p_iETy);

for i=1:length(p_dev_er_rEx)
    p_rel_er_rEx(i)=100*p_dev_er_rEx(i)/abs(p_rETx(i));
    p_rel_er_iEx(i)=100*p_dev_er_iEx(i)/abs(p_iETx(i));
end

for i=1:length(p_dev_er_rEy)
    p_rel_er_rEy(i)=100*p_dev_er_rEy(i)/abs(p_rETy(i));
    p_rel_er_iEy(i)=100*p_dev_er_iEy(i)/abs(p_iETy(i));
end

end