function varargout = elegant_plotter(varargin)
% ELEGANT_PLOTTER MATLAB code for elegant_plotter.fig
%      ELEGANT_PLOTTER, by itself, creates a new ELEGANT_PLOTTER or raises the existing
%      singleton*.
%
%      H = ELEGANT_PLOTTER returns the handle to a new ELEGANT_PLOTTER or the handle to
%      the existing singleton*.
%
%      ELEGANT_PLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ELEGANT_PLOTTER.M with the given input arguments.
%
%      ELEGANT_PLOTTER('Property','Value',...) creates a new ELEGANT_PLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before elegant_plotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to elegant_plotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help elegant_plotter

% Last Modified by GUIDE v2.5 18-Apr-2016 22:45:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @elegant_plotter_OpeningFcn, ...
                   'gui_OutputFcn',  @elegant_plotter_OutputFcn, ...
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


% --- Executes just before elegant_plotter is made visible.
function elegant_plotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to elegant_plotter (see VARARGIN)

% Choose default command line output for elegant_plotter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes elegant_plotter wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(hObject,'Toolbar','figure');

% --- Outputs from this function are returned to the command line.
function varargout = elegant_plotter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

global s_ele sig_x_ele sig_y_ele sig_s_ele sdelta_s_ele enx_s_ele eny_s_ele
global s_ele_Twiss alphax_s_ele alphay_s_ele betax_s_ele betay_s_ele psix_s_ele psiy_s_ele etax_s_ele etay_s_ele
global s_m alphax_m alphay_m betax_m betay_m psix_m psiy_m etax_m etay_m

filename01='beam_sigma_mat_functions.o';
% foramt (s,sigma_x,sigma_y,sigma_s,sigma_delta,emit_norm_x,emit_norm_y)
delimiterIn=' '; headerlinesIn=0;
beam_sigma_mat=importdata(filename01,delimiterIn,headerlinesIn);    
s_ele=beam_sigma_mat(:,1);
sig_x_ele=beam_sigma_mat(:,2);                  % sigma_x in m
sig_y_ele=beam_sigma_mat(:,3);                  % sigma_y in m
sig_s_ele=beam_sigma_mat(:,4);                  % sigma_s in m
sdelta_s_ele=beam_sigma_mat(:,5);
enx_s_ele=beam_sigma_mat(:,6);                  % enx_s in m
eny_s_ele=beam_sigma_mat(:,7);                  % eny_s in m


filename02='Twiss_lattice.o';
% foramt (s,alphax,alphay,betax,betay,psix,psiy)
Twiss_mat=importdata(filename02,delimiterIn,headerlinesIn);    
s_ele_Twiss=Twiss_mat(:,1);
alphax_s_ele=Twiss_mat(:,2);
alphay_s_ele=Twiss_mat(:,3);
betax_s_ele=Twiss_mat(:,4);                     % in m
betay_s_ele=Twiss_mat(:,5);                     % in m
psix_s_ele=Twiss_mat(:,6);                      % in rad
psiy_s_ele=Twiss_mat(:,7);                      % in rad
etax_s_ele=Twiss_mat(:,8);                      % in m
etay_s_ele=Twiss_mat(:,9);                      % in m

%[s_m,alphax_m,alphay_m,betax_m,betay_m,psix_m,psiy_m,etax_m,etay_m]=dipole_info_Twiss(1);
[s_m,alphax_m,alphay_m,betax_m,betay_m,psix_m,psiy_m,etax_m,etay_m]=dipole_info_Twiss_new(1);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('mag_layout.o');

global s_ele sig_x_ele

scale_factor=0.1*(abs(max(sig_x_ele)-min(sig_x_ele)))*1e3;
set(gca,'FontSize',40,'linewidth',5); plot(s_ele,1e3*sig_x_ele,'r','linewidth',3); xlabel('s (m)'); ylabel('\sigma_x(s) (mm)'); hold on; axis('tight');
if (scale_factor>0)
    set(gca,'FontSize',40,'linewidth',5); 
    area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'black', 'EdgeColor', 'black'); hold off; axis('tight'); 
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla;


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('mag_layout.o');

global s_ele sig_y_ele

scale_factor=0.1*(abs(max(sig_y_ele)-min(sig_y_ele)))*1e3;
set(gca,'FontSize',40,'linewidth',5); plot(s_ele,1e3*sig_y_ele,'b','linewidth',3); xlabel('s (m)'); ylabel('\sigma_y(s) (mm)'); hold on; axis('tight');
if (scale_factor>0) 
    set(gca,'FontSize',40,'linewidth',5); 
    area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'black', 'EdgeColor', 'black'); hold off; axis('tight'); 
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('mag_layout.o');

global s_ele sig_s_ele

scale_factor=0.1*(abs(max(sig_s_ele)-min(sig_s_ele)))*1e3;
set(gca,'FontSize',40,'linewidth',5); plot(s_ele,1e3*sig_s_ele,'m','linewidth',3); xlabel('s (m)'); ylabel('\sigma_s(s) (mm)'); hold on; axis('tight');
if (scale_factor>0) 
    set(gca,'FontSize',40,'linewidth',5); 
    area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'black', 'EdgeColor', 'black'); hold off; axis('tight'); 
end

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('mag_layout.o');

global s_ele enx_s_ele

scale_factor=0.1*(abs(max(enx_s_ele)-min(enx_s_ele)))*5e6;
set(gca,'FontSize',40,'linewidth',5); plot(s_ele,1e6*enx_s_ele,'c','linewidth',3); xlabel('s (m)'); ylabel('\epsilon_x(s) (\mum)'); hold on; axis('tight');
if (scale_factor>0)
    set(gca,'FontSize',40,'linewidth',5); 
    area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'black', 'EdgeColor', 'black'); hold off; axis('tight'); 
end

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('mag_layout.o');

global s_ele eny_s_ele

scale_factor=0.1*(abs(max(eny_s_ele)-min(eny_s_ele)))*5e6;
set(gca,'FontSize',40,'linewidth',5); plot(s_ele,1e6*eny_s_ele,'b','linewidth',3); xlabel('s (m)'); ylabel('\epsilon_y(s) (\mum)'); hold on; axis('tight');
if (scale_factor>0)
    set(gca,'FontSize',40,'linewidth',5); 
    area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'black', 'EdgeColor', 'black'); hold off; axis('tight'); 
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('mag_layout.o');

global s_ele_Twiss alphax_s_ele
global s_m alphax_m

scale_factor=0.1*(abs(max(alphax_s_ele)-min(alphax_s_ele)));
set(gca,'FontSize',40,'linewidth',5); plot(s_ele_Twiss,alphax_s_ele,'r','linewidth',3); xlabel('s (m)'); ylabel('\alpha_x(s)'); hold on; axis('tight');
set(gca,'FontSize',40,'linewidth',5); plot(s_m,alphax_m,'bo','linewidth',5); xlabel('s (m)'); ylabel('\alpha_x(s)'); hold on; axis('tight');
if (scale_factor>0) 
    set(gca,'FontSize',40,'linewidth',5); 
    area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'black', 'EdgeColor', 'black'); hold off; axis('tight'); 
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('mag_layout.o');

global s_ele_Twiss alphay_s_ele
global s_m alphay_m

scale_factor=0.1*(abs(max(alphay_s_ele)-min(alphay_s_ele)));
set(gca,'FontSize',40,'linewidth',5); plot(s_ele_Twiss,alphay_s_ele,'g','linewidth',3); xlabel('s (m)'); ylabel('\alpha_y(s)'); hold on; axis('tight');
set(gca,'FontSize',40,'linewidth',5); plot(s_m,alphay_m,'ro','linewidth',5); xlabel('s (m)'); ylabel('\alpha_y(s)'); hold on; axis('tight');
if (scale_factor>0)
    set(gca,'FontSize',40,'linewidth',5); 
    area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'black', 'EdgeColor', 'black'); hold off; axis('tight'); 
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('mag_layout.o');

global s_ele_Twiss betax_s_ele
global s_m betax_m

scale_factor=0.1*(abs(max(betax_s_ele)-min(betax_s_ele)));
set(gca,'FontSize',40,'linewidth',5); plot(s_ele_Twiss,betax_s_ele,'b','linewidth',3); xlabel('s (m)'); ylabel('\beta_x(s) (m)'); hold on; axis('tight');
set(gca,'FontSize',40,'linewidth',5); plot(s_m,betax_m,'ro','linewidth',5); xlabel('s (m)'); ylabel('\beta_x(s) (m)'); hold on; axis('tight');

if (scale_factor>0)
    set(gca,'FontSize',40,'linewidth',5); 
    area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'black', 'EdgeColor', 'black'); hold off; axis('tight'); 
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('mag_layout.o');

global s_ele_Twiss betay_s_ele
global s_m betay_m

scale_factor=0.1*(abs(max(betay_s_ele)-min(betay_s_ele)));
set(gca,'FontSize',40,'linewidth',5); plot(s_ele_Twiss,betay_s_ele,'m','linewidth',3); xlabel('s (m)'); ylabel('\beta_x(s) (m)'); hold on; axis('tight');
set(gca,'FontSize',40,'linewidth',5); plot(s_m,betay_m,'ro','linewidth',5); xlabel('s (m)'); ylabel('\beta_y(s) (m)'); hold on; axis('tight');

if (scale_factor>0)
    set(gca,'FontSize',40,'linewidth',5); 
    area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'black', 'EdgeColor', 'black'); hold off; axis('tight'); 
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('mag_layout.o');

global s_ele_Twiss psix_s_ele
global s_m psix_m

scale_factor=0.1*(abs(max(psix_s_ele)-min(psix_s_ele)));
set(gca,'FontSize',40,'linewidth',5); plot(s_ele_Twiss,psix_s_ele/pi,'c','linewidth',3); xlabel('s (m)'); ylabel('\psi_x(s) (unit of \pi)'); hold on; axis('tight');
set(gca,'FontSize',40,'linewidth',5); plot(s_m,psix_m/pi,'ro','linewidth',5); xlabel('s (m)'); ylabel('\psi_x(s) (unit of \pi)'); hold on; axis('tight');

if (scale_factor>0)
    set(gca,'FontSize',40,'linewidth',5); 
    area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'black', 'EdgeColor', 'black'); hold off; axis('tight'); 
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('mag_layout.o');

global s_ele_Twiss psiy_s_ele
global s_m psiy_m

scale_factor=0.1*(abs(max(psiy_s_ele)-min(psiy_s_ele)));
set(gca,'FontSize',40,'linewidth',5); plot(s_ele_Twiss,psiy_s_ele/pi,'b','linewidth',3); xlabel('s (m)'); ylabel('\psi_y(s) (unit of \pi)'); hold on; axis('tight');
set(gca,'FontSize',40,'linewidth',5); plot(s_m,psiy_m/pi,'ro','linewidth',5); xlabel('s (m)'); ylabel('\psi_y(s) (unit of \pi)'); hold on; axis('tight');

if (scale_factor>0)
    set(gca,'FontSize',40,'linewidth',5); 
    area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'black', 'EdgeColor', 'black'); hold off; axis('tight'); 
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('mag_layout.o');

global s_ele sdelta_s_ele

scale_factor=0.1*(abs(max(sdelta_s_ele)-min(sdelta_s_ele)));

set(gca,'FontSize',40,'linewidth',5); plot(s_ele,sdelta_s_ele,'r','linewidth',3); xlabel('s (m)'); ylabel('\sigma_{\delta}(s)'); hold on; axis('tight');
if (scale_factor>0)
    set(gca,'FontSize',40,'linewidth',5); 
    area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'black', 'EdgeColor', 'black'); hold off; axis('tight');
end


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('mag_layout.o');

global s_ele_Twiss etax_s_ele
global s_m etax_m

scale_factor=0.1*(abs(max(etax_s_ele)-min(etax_s_ele)));

set(gca,'FontSize',40,'linewidth',5); plot(s_ele_Twiss,etax_s_ele,'r','linewidth',3); xlabel('s (m)'); ylabel('\eta_x(s) (m)'); hold on; axis('tight');
set(gca,'FontSize',40,'linewidth',5); plot(s_m,etax_m,'ro','linewidth',5); xlabel('s (m)'); ylabel('\eta_x(s) (m)'); hold on; axis('tight');

if (scale_factor>0)
    set(gca,'FontSize',40,'linewidth',5); 
    area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'black', 'EdgeColor', 'black'); hold off; axis('tight');
end

% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('mag_layout.o');

global s_ele_Twiss etay_s_ele
global s_m etay_m

scale_factor=0.1*(abs(max(etay_s_ele)-min(etay_s_ele)));

set(gca,'FontSize',40,'linewidth',5); plot(s_ele_Twiss,etay_s_ele,'r','linewidth',3); xlabel('s (m)'); ylabel('\eta_y(s) (m)'); hold on; axis('tight');
set(gca,'FontSize',40,'linewidth',5); plot(s_m,etay_m,'ro','linewidth',5); xlabel('s (m)'); ylabel('\eta_y(s) (m)'); hold on; axis('tight');

if (scale_factor>0)
    set(gca,'FontSize',40,'linewidth',5); 
    area(mag_layout(:,1),scale_factor*mag_layout(:,2),'FaceColor', 'black', 'EdgeColor', 'black'); hold off; axis('tight');
end
