function f_plotSurface(nodeMear,fname)
nodeCount = size(nodeMear,1);
for i=1:nodeCount
    eval(['!wb_command -gifti-label-to-roi surfFS.lh.MBM_cortex_vM.label.gii -key ' num2str(i) ' L' num2str(i) '.func.gii']);
    if i==1
        g_temp = gifti(['L' num2str(i) '.func.gii']);
        g_temp.cdata(g_temp.cdata==1)=nodeMear(i);
    end
        g_temp_next = gifti(['L' num2str(i) '.func.gii']);
        g_temp_next.cdata(g_temp_next.cdata==1)=nodeMear(i);
        g_temp.cdata = g_temp.cdata + g_temp_next.cdata;
    if i==nodeCount
        save(g_temp,[fname,'.func.gii'],'Base64Binary');
    end
end
!rm -rf L*.func.gii
end

