

function varargout = volterra_plotter(varargin)
% VOLTERRA_PLOTTER MATLAB code for volterra_plotter.fig
%      VOLTERRA_PLOTTER, by itself, creates a new VOLTERRA_PLOTTER or raises the existing
%      singleton*.
%
%      H = VOLTERRA_PLOTTER returns the handle to a new VOLTERRA_PLOTTER or the handle to
%      the existing singleton*.
%
%      VOLTERRA_PLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VOLTERRA_PLOTTER.M with the given input arguments.
%
%      VOLTERRA_PLOTTER('Property','Value',...) creates a new VOLTERRA_PLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before volterra_plotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to volterra_plotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help volterra_plotter

% Last Modified by GUIDE v2.5 16-May-2015 20:42:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @volterra_plotter_OpeningFcn, ...
                   'gui_OutputFcn',  @volterra_plotter_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before volterra_plotter is made visible.
function volterra_plotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to volterra_plotter (see VARARGIN)

% Choose default command line output for volterra_plotter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes volterra_plotter wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(hObject,'Toolbar','figure');

%{
load('mag_layout.o'); warning('off');
global mag_layout;

load('workplace.mat'); warning('off');
global s G G_zero G_first G_second G_third G_fourth G_fifth G_sixth G_seventh G_eighth G_ninth G_tenth G_c_matrix Gfp Gp;
global Gf Gf_zero Gf_first Gf_second Gf_third Gf_fourth Gf_fifth Gf_sixth Gf_seventh Gf_eighth Gf_ninth Gf_tenth;
global I_b lambda_array scan_num iiterative ilast find_TWLA egamma_vec betax0 betay0 start_pos end_pos;
global s_ele R16_ele R36_ele R51_ele R52_ele R53_ele R54_ele R56_ele C_ele;
%}

% --- Outputs from this function are returned to the command line.
function varargout = volterra_plotter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global s G G_zero G_first G_second G_third G_fourth G_fifth G_sixth G_seventh G_eighth G_ninth G_tenth;
%global scan_num iiterative ilast;
load('workplace.mat');
if ((scan_num==1))

    
    if (iiterative==1 && ilast==0)
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_zero,'r-','linewidth',3);    xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_first,'y-','linewidth',3);   xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_second,'g-','linewidth',3);  xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_third,'c-','linewidth',3);   xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_fourth,'b-','linewidth',3);  xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_fifth,'m-','linewidth',3);   xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_sixth,'r--','linewidth',3);  xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_seventh,'y--','linewidth',3);xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_eighth,'g--','linewidth',3); xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_ninth,'c--','linewidth',3);  xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_tenth,'b--','linewidth',3);  xlabel('s (m)'); ylabel('G(s)'); hold off;grid off; axis('tight');
    else
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G,'k-','linewidth',3);         xlabel('s (m)'); ylabel('G(s)');  grid off; axis('tight');
    end
else
    msgbox('scan_num > 1, can only plot gain spectrum or gain map...','','warn');
end



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global lambda_array Gf Gf_zero Gf_first Gf_second Gf_third Gf_fourth Gf_fifth Gf_sixth Gf_seventh Gf_eighth Gf_ninth Gf_tenth;
%global scan_num iiterative;
load('workplace.mat');
if (scan_num>=2)
        if (iiterative==1)
            set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,Gf,'k-','linewidth',3);         xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold on; grid off; axis('tight');
            set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,Gf_zero,'r-','linewidth',3);    xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold on; grid off; axis('tight');
            set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,Gf_first,'y-','linewidth',3);   xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold on; grid off; axis('tight');
            set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,Gf_second,'g-','linewidth',3);  xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold on; grid off; axis('tight');
            set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,Gf_third,'c-','linewidth',3);   xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold on; grid off; axis('tight');
            set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,Gf_fourth,'b-','linewidth',3);  xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold on; grid off; axis('tight');
            set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,Gf_fifth,'m-','linewidth',3);   xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold on; grid off; axis('tight');
            set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,Gf_sixth,'r--','linewidth',3);  xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold on; grid off; axis('tight');
            set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,Gf_seventh,'y--','linewidth',3);xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold on; grid off; axis('tight');
            set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,Gf_eighth,'g--','linewidth',3); xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold on; grid off; axis('tight');
            set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,Gf_ninth,'c--','linewidth',3);  xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold on; grid off; axis('tight');
            set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,Gf_tenth,'b--','linewidth',3);  xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold off;grid off; axis('tight');
        else
            set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,Gf,'k-','linewidth',3); xlabel('\lambda (\mum)'); ylabel('G_{f}');  grid off; axis('tight');
        end
else
    msgbox('scan_num = 1, can only plot gain function...','','warn');
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
%global s_ele R16_ele R36_ele R51_ele R52_ele R53_ele R54_ele R56_ele;
%global find_TWLA egamma_vec;
%global mag_layout;
load('workplace.mat'); load('mag_layout.o');

a=get(gcbo,'Value');
if (a==1)
    scale_factor=0.1*(abs(max(R16_ele)-min(R16_ele)));
    set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R16_ele,'b-','linewidth',3); xlabel('s (m)'); ylabel('R_{16} (cm)'); hold on;     grid off;  axis('tight');
    set(gca,'FontSize',40,'linewidth',5); area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'r', 'EdgeColor', 'r'); hold off; axis('tight');
end
if (a==2)
    scale_factor=0.1*(abs(max(R36_ele)-min(R36_ele)));
    set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R36_ele,'b-','linewidth',3); xlabel('s (m)'); ylabel('R_{36} (cm)');  hold on;    grid off;  axis('tight');
    set(gca,'FontSize',40,'linewidth',5); area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'r', 'EdgeColor', 'r'); hold off; axis('tight');
end
if (a==3)
    scale_factor=0.1*(abs(max(R51_ele)-min(R51_ele)));
    set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R51_ele,'b-','linewidth',3); xlabel('s (m)'); ylabel('R_{51}');     hold on;      grid off;  axis('tight');
    set(gca,'FontSize',40,'linewidth',5); area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'r', 'EdgeColor', 'r'); hold off; axis('tight');
end
if (a==4)
    scale_factor=0.1*(abs(max(R52_ele)-min(R52_ele)));
    set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R52_ele,'b-','linewidth',3); xlabel('s (m)'); ylabel('R_{52} (cm)');  hold on;    grid off;  axis('tight');
    set(gca,'FontSize',40,'linewidth',5); area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'r', 'EdgeColor', 'r'); hold off; axis('tight');
end
if (a==5)
    scale_factor=0.1*(abs(max(R53_ele)-min(R53_ele)));
    %if (scale_factor==0); scale_factor=0.1; end;
    set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R53_ele,'b-','linewidth',3); xlabel('s (m)'); ylabel('R_{53}');     hold on;      grid off;  axis('tight');
    set(gca,'FontSize',40,'linewidth',5); area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'r', 'EdgeColor', 'r'); hold off; axis('tight');
end
if (a==6)
    scale_factor=0.1*(abs(max(R54_ele)-min(R54_ele)));
    %if (scale_factor==0); scale_factor=0.1; end;
    set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R54_ele,'b-','linewidth',3); xlabel('s (m)'); ylabel('R_{54} (cm)');  hold on;    grid off;  axis('tight');
    set(gca,'FontSize',40,'linewidth',5); area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'r', 'EdgeColor', 'r'); hold off; axis('tight');
end
if (a==7)
    if (find_TWLA==1)
        scale_factor=0.1*(abs(max(R56_ele.*egamma_vec)-min(R56_ele.*egamma_vec)));
    	set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R56_ele.*egamma_vec,'b-','linewidth',3); xlabel('s (m)'); ylabel('R_{56} (cm)');  hold on;    grid off;  axis('tight');
        set(gca,'FontSize',40,'linewidth',5); area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'r', 'EdgeColor', 'r'); hold off; axis('tight');
    else
        scale_factor=0.1*(abs(max(R56_ele)-min(R56_ele)));
        set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R56_ele,'b-','linewidth',3); xlabel('s (m)'); ylabel('R_{56} (cm)');   hold on;   grid off;  axis('tight');
        set(gca,'FontSize',40,'linewidth',5); area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'r', 'EdgeColor', 'r'); hold off; axis('tight');
    end
end



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global s_ele C_ele;
%global mag_layout;
load('workplace.mat'); load('mag_layout.o');
scale_factor=0.1*(abs(max(C_ele)-min(C_ele)));
if (scale_factor==0); scale_factor=0.5; end;
set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,C_ele,'b-','linewidth',3); xlabel('s (m)'); ylabel('Compression factor'); grid off; hold on; axis('tight');
set(gca,'FontSize',40,'linewidth',5); area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'r', 'EdgeColor', 'r'); hold off; axis('tight');

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global s_ele C_ele I_b;
%global mag_layout;
load('workplace.mat'); load('mag_layout.o');
scale_factor=0.1*(abs(max(C_ele*I_b)-min(C_ele*I_b)));
if (scale_factor==0); scale_factor=0.5; end;
set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,C_ele*I_b,'b-','linewidth',3); xlabel('s (m)'); ylabel('Peak current (A)'); grid off; hold on; axis('tight');
set(gca,'FontSize',40,'linewidth',5); area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'r', 'EdgeColor', 'r'); hold off; axis('tight');


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global s_ele R16_ele R36_ele R51_ele R52_ele R53_ele R54_ele R56_ele C_ele betax0 betay0 start_pos end_pos;
load('workplace.mat');
msgbox('This takes a second...','','warn');
quilt_plot;

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global s scan_num lambda_array G_c_matrix;
load('workplace.mat');
if (scan_num>=2)
	figure(123); 
    set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,abs(G_c_matrix)); view(-38,32); xlabel('\lambda (\mum)','FontSize',40); ylabel('s (m)','FontSize',40); zlabel('G(s,\lambda)','FontSize',40); axis('tight'); 
    shading('interp'); rotate3d on; h=colorbar(); set(h,'FontSize',30);
else
    msgbox('scan_num = 1, can only plot gain function...','','warn');
end

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global scan_num lambda_array Gfp;
load('workplace.mat');
if (scan_num>=2)
	set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,Gfp,'k-','linewidth',3); xlabel('\lambda (\mum)'); ylabel('G_{f}^{p}');  grid off; axis('tight');
	%[hAx,hLine1,hLine2]=plotyy(lambda_array*10^4,Gf,lambda_array*10^4,Gfp); grid off; set(hAx,'FontSize',40,'linewidth',5); set(hLine1,'linewidth',5); set(hLine2,'linewidth',5); hLine1.LineStyle='r-'; hLine2.LineStyle='b-'; hLine1.Linewidth=5; hLine2.Linewidth=5; xlabel('\lambda (\mum)'); ylabel(hAx(1),'G_{f}'); ylabel(hAx(2),'G_{f}^{p}');
    msgbox('relative quantity, under construction...','','warn');
else
    msgbox('scan_num = 1, can only plot energy modulation function...','','warn');
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global scan_num s Gp;
load('workplace.mat');
if ((scan_num==1))
    set(gca,'FontSize',40,'linewidth',5); plot(s/100,Gp,'k-','linewidth',3);    xlabel('s (m)'); ylabel('G^{p}(s)');  grid off; axis('tight'); hold off;
	%[hAx,hLine1,hLine2]=plotyy(s/100,G,s/100,Gp); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40,'linewidth',5); set(hLine1,'linewidth',5); set(hLine2,'linewidth',5); hLine1.LineStyle='r-'; hLine2.LineStyle='b-'; hLine1.Linewidth=5; hLine2.Linewidth=5; xlabel('s (m)'); ylabel(hAx(1),'G(s)'); ylabel(hAx(2),'G^{p}(s)');
    msgbox('relative quantity, under construction...','','warn');
else
    msgbox('scan_num > 1, can only plot energy modulation spectrum...','','warn');
end        


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global scan_num s Gp;
%global mag_layout;
load('workplace.mat'); load('mag_layout.o');
scale_factor=0.1*(abs(max(Gp)-min(Gp)));
if ((scan_num==1))
    set(gca,'FontSize',40,'linewidth',5); plot(s/100,Gp,'k-','linewidth',3);    xlabel('s (m)'); ylabel('G^{p}(s)'); hold on; grid off; axis('tight');
    set(gca,'FontSize',40,'linewidth',5); area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'r', 'EdgeColor', 'r'); axis('tight'); hold off;
    msgbox('relative quantity, under construction...','','warn');
else
    msgbox('scan_num > 1, can only plot energy modulation spectrum...','','warn');
end  

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global s G G_zero G_first G_second G_third G_fourth G_fifth G_sixth G_seventh G_eighth G_ninth G_tenth;
%global scan_num iiterative ilast;
%global mag_layout;
load('workplace.mat'); load('mag_layout.o');
scale_factor=0.1*(abs(max(G)-min(G)));
if ((scan_num==1))

    if (iiterative==1 && ilast==0)
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_zero,'r-','linewidth',3);    xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_first,'y-','linewidth',3);   xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_second,'g-','linewidth',3);  xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_third,'c-','linewidth',3);   xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_fourth,'b-','linewidth',3);  xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_fifth,'m-','linewidth',3);   xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_sixth,'r--','linewidth',3);  xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_seventh,'y--','linewidth',3);xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_eighth,'g--','linewidth',3); xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_ninth,'c--','linewidth',3);  xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_tenth,'b--','linewidth',3);  xlabel('s (m)'); ylabel('G(s)'); hold on ;grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G,'k-','linewidth',3);         xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'r', 'EdgeColor', 'r'); hold off; axis('tight');
    else
        set(gca,'FontSize',40,'linewidth',5); plot(s/100,G,'k-','linewidth',3);         xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        set(gca,'FontSize',40,'linewidth',5); area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'r', 'EdgeColor', 'r'); hold off; axis('tight');
    end
else
    msgbox('scan_num > 1, can only plot gain spectrum or gain map...','','warn');
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
axes(hObject);


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2
%axes(hObject);
imshow('vt-logo.gif');


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes4
%axes(hObject);
imshow('JLab-logo.jpg');

