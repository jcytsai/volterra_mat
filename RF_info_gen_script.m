%clear all;
format long;

global find_TWLA;

%============ load the ELEGANT element geometrical information ===========%
filename1='elements_geom.o';
% format (ElementType,s,theta,phi)
delimiterIn='\t'; headerlinesIn=0;
RF_info=importdata(filename1,delimiterIn,headerlinesIn);

[r,~]=size(RF_info.data);
%=========================================================================%


%============ select elements and read cavity model information ==========%
n=0; p=0; q=0; w=0;
for m=1:1:r
    
    if ((strcmp(RF_info.textdata(m),'RFCA')==1)||...
        (strcmp(RF_info.textdata(m),'RFCW')==1)||...
        (strcmp(RF_info.textdata(m),'TWLA')==1))
        w=w+1;    
        if (w==1)
            n=n+1;
            RF_ent_row_index(n)=m-1;
        end
    else
        if (n~=0 && w~=0)
            p=p+1;
            RF_ext_row_index(p)=m-1;
            w=0;
        end
    end
    
end
%=========================================================================%
[rr,cc]=size(RF_ent_row_index); n=cc;

RF_ent_ext_row_index(1:2:2*cc-1)=RF_ent_row_index;
RF_ent_ext_row_index(2:2:2*cc)=RF_ext_row_index;

%----------------------%
[rrr,ccc]=size(RF_ent_ext_row_index);
if (ccc ~= 2*cc)
    fprintf('warning: something wrong in RF_elements_info.o...\n');
    pause;
end
%----------------------%

RF_ent_s=RF_info.data(RF_ent_row_index,1);
RF_exit_s=RF_info.data(RF_ext_row_index,1);

%============ format conversion ==========================================%
RF_s(1:2:2*n-1)=RF_ent_s;
RF_s(2:2:2*n)=RF_exit_s;
RF_s=RF_s';
%=========================================================================%

%============ load the ELEGANT element parametrical information ==========%
filename2='elements_param.o';
% format (ElementName,ElementType,ParameterName,ParameterValue)
delimiterIn=' '; headerlinesIn=0;
RF_para=importdata(filename2,delimiterIn,headerlinesIn);

[r,~]=size(RF_para.data);
%=========================================================================%


%===== select elements, properties and read cavity model parameters ======%
n=0; p=0; q=0; w=0; m=1;

if (find_TWLA==2)
while (m<=r) % RFCA-type
    if ((strcmp(RF_para.textdata(m,2),'RFCA')==1)||...
        (strcmp(RF_para.textdata(m,2),'RFCW')==1))
            n=n+1;
            RF_data_info(n,1)=RF_para.data(m);               % L (m)
            RF_data_info(n,2)=RF_para.data(m+1)/1e6;         % VOLT (MV)
            RF_data_info(n,3)=RF_para.data(m+2)-90;          % PHASE (deg), use cosine convention. Note: ELEGANT use sine convention.
            RF_data_info(n,4)=RF_para.data(m+3);             % FREQ (Hz)
            m=m+18;
    else
        m=m+1;
    end    
end
else
while (m<=r) % TWLA-type
    if ((strcmp(RF_para.textdata(m,2),'TWLA')==1))
            n=n+1;
            RF_data_info(n,1)=RF_para.data(m);               % L (m)
            RF_data_info(n,2)=RF_para.data(m+4)/1e6;         % EZ (MV/m)
            RF_data_info(n,3)=RF_para.data(m+2)/pi*180*(-1); % PHASE (deg)
            RF_data_info(n,4)=RF_para.data(m+1);             % FREQ (Hz)
            m=m+18;
    else
        m=m+1;
    end    
end
end

[r,~]=size(RF_data_info);
for m=1:1:r
    RF_ddata_info(2*m-1,:)=RF_data_info(m,:);
    RF_ddata_info(2*m,:)=RF_data_info(m,:);
end

[r1,~]=size(RF_s);
[r2,~]=size(RF_ddata_info);

if (r1~=r2)
    fprintf('warning: RF information may be insufficient...\n');
end
%=========================================================================%
RF_mat=[RF_s RF_ddata_info];
%=========================================================================%

%============ summary and output corrected lattice file ==================%

fprintf('generating RF element info ...'); fprintf('done... \n');
fprintf('total number of %d RF elements (RFCA, RFCW and/or TWLA) detected... \n',r);

fileID1 = fopen('RF_elements_info.o','w');
fprintf(fileID1,'%5.5f %5.5f %5.5f %5.5f %5.5f \n',RF_mat');
fclose(fileID1);

%=========================================================================%







