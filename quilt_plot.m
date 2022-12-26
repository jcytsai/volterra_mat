function quilt_plot

global s_ele R16_ele R36_ele R51_ele R52_ele R53_ele R54_ele R56_ele C_ele betax0 betay0 start_pos end_pos;

% this program generates R51(s,s'), R52(s,s'), and R56(s,s') composite 
% tranport function according to the following relations:
% R51(s,s') = C(s)R51(s) - C(s')R51(s')
% R52(s,s') = C(s)R52(s) - C(s')R52(s')
% R56(s,s') = C(s)R56(s) - C(s')R56(s')
% and R56(s'->s) = R56(s) - R56(s') + R51(s')*R52(s) - R51(s)*R52(s')
% where s' is the emission/source point and s is observation/test point
% output results are set in MKS units

% read lattice_transport_functions.o file
% (s,R16,R36,R51,R52,R53,R54,R56)

filename='lattice_transport_functions.o';
mesh_num=1000;

s=linspace(start_pos,end_pos,mesh_num);
tic;

for p=1:1:mesh_num                                       % index for s'
    s_source(p)=s(p);
    %fprintf('%d source step done... \n',p);
    %for q=1:1:mesh_num                                  % index for s
        %s_test=s(q);
        s_test=s;
        tmp01=interp1(s_ele,R51_ele,s_test);             %R51(s)
        tmp02=interp1(s_ele,R51_ele,s_source(p));        %R51(s')
        
        tmp03=interp1(s_ele,R52_ele,s_test);             %R52(s)
        tmp04=interp1(s_ele,R52_ele,s_source(p));        %R52(s')
        %{
        if (ilattice==1)
            tmp11=interp1(s_ele,R53_ele,s_test);         %R53(s)
            tmp12=interp1(s_ele,R53_ele,s_source(p));    %R53(s')
        
            tmp13=interp1(s_ele,R54_ele,s_test);         %R54(s)
            tmp14=interp1(s_ele,R54_ele,s_source(p));    %R54(s')
        end
        %}
        tmp05=interp1(s_ele,R56_ele,s_test);             %R56(s)
        tmp06=interp1(s_ele,R56_ele,s_source(p));        %R56(s')
        
        tmp07=interp1(s_ele,C_ele,s_test);               %C(s)
        tmp08=interp1(s_ele,C_ele,s_source(p));          %C(s')
        %{
        if (ilattice==1)
            R51_mod(p,:)=tmp07.*tmp01-tmp08*tmp02;
            R52_mod(p,:)=tmp07.*tmp03-tmp08*tmp04;
            R53_mod(p,:)=tmp07.*tmp11-tmp08*tmp12;
            R54_mod(p,:)=tmp07.*tmp13-tmp08*tmp14;
            R56_mod(p,:)=tmp07.*tmp05-tmp08*tmp06;
        
            R51_mod(p,1:p)=zeros(1,p);
            R52_mod(p,1:p)=zeros(1,p);
            R53_mod(p,1:p)=zeros(1,p);
            R54_mod(p,1:p)=zeros(1,p);
            R56_mod(p,1:p)=zeros(1,p);
        else
        %}
            R51_mod(p,:)=tmp07.*tmp01-tmp08*tmp02;
            R52_mod(p,:)=tmp07.*tmp03-tmp08*tmp04;
            R56_mod(p,:)=tmp07.*tmp05-tmp08*tmp06;
        
            R51_mod(p,1:p)=zeros(1,p);
            R52_mod(p,1:p)=zeros(1,p);
            R56_mod(p,1:p)=zeros(1,p);
        %end
        
        %tmp05 => tmp01=interp1(s_ele,R56_ele,s_test);
        %tmp06 => tmp02=interp1(s_ele,R56_ele,s_source(p));
        %tmp02*tmp03 => tmp03=interp1(s_ele,R51_ele,s_source(p))*interp1(s_ele,R52_ele,s_test);
        %tmp01*tmp04 => tmp04=interp1(s_ele,R51_ele,s_test)*interp1(s_ele,R52_ele,s_source(p));
        %{
        if (ilattice==1)
            tmp15=tmp12*tmp13;
            tmp16=tmp11*tmp14;
            R56_map(p,:)=tmp05-tmp06+tmp02*tmp03-tmp01*tmp04+tmp15-tmp16;
            R56_map(p,1:p)=zeros(1,p);
        else
        %}
            R56_map(p,:)=tmp05-tmp06+tmp02*tmp03-tmp01*tmp04;
            R56_map(p,1:p)=zeros(1,p);
        %end
        
    %end
end
toc;

figure(99);
%{
    set(gca,'FontSize',40); subplot(3,2,1); h1=mesh(s_test/100,s_source/100,betax0/100*R51_mod.^2); view(0,90); grid off; axis('tight'); ylabel('s\prime source (m)','FontSize',40); colorbar;
    set(gca,'FontSize',40); subplot(3,2,2); h2=mesh(s_test/100,s_source/100,(1/(betax0/100))*(R52_mod/100).^2); grid off; view(0,90); axis('tight'); colorbar;
    set(gca,'FontSize',40); subplot(3,2,3); h=mesh(s_test/100,s_source/100,betay0/100*R53_mod.^2); view(0,90); grid off; axis('tight'); ylabel('s\prime source (m)','FontSize',40); colorbar;
    set(gca,'FontSize',40); subplot(3,2,4); hh=mesh(s_test/100,s_source/100,(1/(betay0/100))*(R54_mod/100).^2); grid off; view(0,90); axis('tight'); colorbar;
    set(gca,'FontSize',40); subplot(3,2,5); h3=mesh(s_test/100,s_source/100,(R56_mod/100).^2); view(0,90); grid off; axis('tight'); ylabel('s\prime source (m)','FontSize',40); xlabel('s test (m)','FontSize',40); colorbar;
    set(gca,'FontSize',40); subplot(3,2,6); h4=mesh(s_test/100,s_source/100,abs(R56_map)/100); view(0,90); grid off; axis('tight'); xlabel('s test (m)','FontSize',40); colorbar;
%}
    set(gca,'FontSize',40);subplot(2,2,1); h1=mesh(s_test/100,s_source/100,betax0/100*R51_mod.^2);             view(0,90); grid off; axis('tight'); ylabel('s\prime source (m)','FontSize',40); colorbar;
    set(gca,'FontSize',40);subplot(2,2,2); h2=mesh(s_test/100,s_source/100,(1/(betax0/100))*(R52_mod/100).^2); view(0,90); grid off; axis('tight'); colorbar;
    set(gca,'FontSize',40);subplot(2,2,3); h3=mesh(s_test/100,s_source/100,(R56_mod/100).^2);                  view(0,90); grid off; axis('tight'); ylabel('s\prime source (m)','FontSize',40); xlabel('s test (m)','FontSize',40); colorbar;
    set(gca,'FontSize',40);subplot(2,2,4); h4=mesh(s_test/100,s_source/100,abs(R56_map)/100);                  view(0,90); grid off; axis('tight'); set(gca,'FontSize',40); xlabel('s test (m)','FontSize',40); h=colorbar; set(h,'FontSize',40);





