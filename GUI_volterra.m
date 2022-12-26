
function varargout = GUI_volterra(varargin)
% GUI_VOLTERRA MATLAB code for GUI_volterra.fig
%      GUI_VOLTERRA, by itself, creates a new GUI_VOLTERRA or raises the existing
%      singleton*.
%
%      H = GUI_VOLTERRA returns the handle to a new GUI_VOLTERRA or the handle to
%      the existing singleton*.
%
%      GUI_VOLTERRA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_VOLTERRA.M with the given input arguments.
%
%      GUI_VOLTERRA('Property','Value',...) creates a new GUI_VOLTERRA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_volterra_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_volterra_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_volterra

% Last Modified by GUIDE v2.5 11-Oct-2022 13:24:43

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_volterra_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_volterra_OutputFcn, ...
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


% --- Executes just before GUI_volterra is made visible.
function GUI_volterra_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_volterra (see VARARGIN)

% Choose default command line output for GUI_volterra
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_volterra wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_volterra_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

global energy I_b C_factor emit_norm_x emit_norm_y betax0 betay0 alphax0 alphay0 sigma_delta chirp start_pos end_pos lambda_start01 lambda_end01 scan_num01 mesh_num lambda_start02 lambda_end02 scan_num02 iplot_lattice iplot_gain_func iplot_gain_spec isave;
global iiterative ilast iCSR_ss iCSR_tr iCSR_drift iLSC ilinac iGaussian iEnergy_mod_calc iEnergy_gain_calc iplot_energy_mod surfplot_gain_spec_func ssCSR_model LSC_model issCSRpp iquilt_plot idm_analysis iplot_I_b full_pipe_height iRF_ele;
global iIBS DO_X DO_Y DO_Z iDz iDzz itransLD use_orig_mat_output plot_Zk_vs_s iSES CSRdrifmodel Derbenev gen_CSR_LSC_param checkbox31 checkbox32

energy_s=num2str(energy); set(handles.energy,'String',energy_s);
I_b_s=num2str(I_b); set(handles.I_b,'String',I_b_s);
C_factor_s=num2str(C_factor); set(handles.C_factor,'String',C_factor_s);
emit_norm_x_s=num2str(emit_norm_x); set(handles.emit_norm_x,'String',emit_norm_x_s);
emit_norm_y_s=num2str(emit_norm_y); set(handles.emit_norm_y,'String',emit_norm_y_s);
betax0_s=num2str(betax0); set(handles.betax0,'String',betax0_s);
betay0_s=num2str(betay0); set(handles.betay0,'String',betay0_s);
alphax0_s=num2str(alphax0); set(handles.alphax0,'String',alphax0_s);
alphay0_s=num2str(alphay0); set(handles.alphay0,'String',alphay0_s);
sigma_delta_s=num2str(sigma_delta); set(handles.sigma_delta,'String',sigma_delta_s);
chirp_s=num2str(chirp); set(handles.chirp,'String',chirp_s);
start_pos_s=num2str(start_pos); set(handles.start_pos,'String',start_pos_s);
end_pos_s=num2str(end_pos); set(handles.end_pos,'String',end_pos_s);
lambda_start01_s=num2str(lambda_start01); set(handles.lambda_start01,'String',lambda_start01_s);
lambda_end01_s=num2str(lambda_end01); set(handles.lambda_end01,'String',lambda_end01_s);
scan_num01_s=num2str(scan_num01); set(handles.scan_num01,'String',scan_num01_s);
lambda_start02_s=num2str(lambda_start02); set(handles.lambda_start02,'String',lambda_start02_s);
lambda_end02_s=num2str(lambda_end02); set(handles.lambda_end02,'String',lambda_end02_s);
scan_num02_s=num2str(scan_num02); set(handles.scan_num02,'String',scan_num02_s);
mesh_num_s=num2str(mesh_num); set(handles.mesh_num,'String',mesh_num_s);
%iplot_lattice_s=num2str(iplot_lattice); set(handles.iplot_lattice,'String',iplot_lattice_s);
%iplot_gain_func_s=num2str(iplot_gain_func); set(handles.iplot_gain_func,'String',iplot_gain_func_s);
%iplot_gain_spec_s=num2str(iplot_gain_spec); set(handles.iplot_gain_spec,'String',iplot_gain_spec_s);
%iplot_energy_mod_s=num2str(iplot_energy_mod); set(handles.iplot_energy_mod,'String',iplot_energy_mod_s);
%isave_s=num2str(isave); set(handles.isave,'String',isave_s);
%iiterative_s=num2str(iiterative); set(handles.iiterative,'String',iiterative_s);
%ilast_s=num2str(ilast); set(handles.ilast,'String',ilast_s);
%iCSR_tr_s=num2str(iCSR_tr); set(handles.iCSR_tr,'String',iCSR_tr_s);
%iCSR_ss_s=num2str(iCSR_ss); set(handles.iCSR_ss,'String',iCSR_ss_s);
%iCSR_drift_s=num2str(iCSR_drift); set(handles.iCSR_drift,'String',iCSR_drift_s);
%iLSC_s=num2str(iLSC); set(handles.iLSC,'String',iLSC_s);
%ilinac_s=num2str(ilinac); set(handles.ilinac,'String',ilinac_s);
%iGaussian_s=num2str(iGaussian); set(handles.iGaussian,'String',iGaussian_s);
%iEnergy_mod_calc_s=num2str(iEnergy_mod_calc); set(handles.iEnergy_mod_calc,'String',iEnergy_mod_calc_s);
%iEnergy_gain_calc_s=num2str(iEnergy_gain_calc); set(handles.iEnergy_gain_calc,'String',iEnergy_gain_calc_s);
%idm_analysis_s=num2str(idm_analysis); set(handles.idm_analysis,'String',idm_analysis_s);
iplot_lattice=handles.iplot_lattice.Value;
iplot_gain_func=handles.iplot_gain_func.Value;
iplot_gain_spec=handles.iplot_gain_spec.Value;
iplot_energy_mod=handles.iplot_energy_mod.Value;
%isave=handles.isave.Value;
iiterative=handles.iiterative.Value;
ilast=handles.ilast.Value;
iCSR_tr=handles.iCSR_tr.Value;
iCSR_ss=handles.iCSR_ss.Value;
iCSR_drift=handles.iCSR_drift.Value;
iLSC=handles.iLSC.Value;
%ilinac=handles.ilinac.Value;
%iGaussian=handles.iGaussian.Value;
iEnergy_mod_calc=handles.iEnergy_mod_calc.Value;
iEnergy_gain_calc=handles.iEnergy_gain_calc.Value;
idm_analysis=handles.idm_analysis.Value;
iquilt_plot=handles.iquilt_plot.Value;
iplot_I_b=handles.iplot_I_b.Value;
%iRF_ele=handles.iRF_ele.Value;
iIBS=handles.iIBS.Value;
DO_X=handles.DO_X.Value;
DO_Y=handles.DO_Y.Value;
DO_Z=handles.DO_Z.Value;
iDz=handles.iDz.Value;
iDzz=handles.iDzz.Value;
itransLD=handles.itransLD.Value;
iSES=handles.iSES.Value;
%use_orig_mat_output=handles.use_orig_mat_output.Value;
plot_Zk_vs_s=handles.plot_Zk_vs_s.Value;
Derbenev=handles.Derbenev.Value;

%surfplot_gain_spec_func_s=num2str(surfplot_gain_spec_func); set(handles.surfplot_gain_spec_func,'String',surfplot_gain_spec_func_s);
%issCSRpp_s=num2str(issCSRpp); set(handles.issCSRpp,'String',issCSRpp_s);
surfplot_gain_spec_func=handles.surfplot_gain_spec_func.Value;
issCSRpp=handles.issCSRpp.Value;

ssCSR_model=handles.popupmenu6.Value;
CSRdrifmodel=handles.popupmenu7.Value;
LSC_model=handles.popupmenu8.Value;
gen_CSR_LSC_param=handles.gen_CSR_LSC_param.Value;
checkbox31=handles.checkbox31.Value;  % generate CSRCSBEND, CSRDRIFT, LSCDRIFT
checkbox32=handles.checkbox32.Value;  % save as SDDS

%ssCSR_model_s=num2str(ssCSR_model); set(handles.ssCSR_model,'String',ssCSR_model_s);
%CSRdrifmodel_s=num2str(CSRdrifmodel); set(handles.CSRdrifmodel,'String',CSRdrifmodel_s);
%LSC_model_s=num2str(LSC_model); set(handles.LSC_model,'String',LSC_model_s);
%iquilt_plot_s=num2str(iquilt_plot); set(handles.iquilt_plot,'String',iquilt_plot_s);
%iplot_I_b_s=num2str(iplot_I_b); set(handles.iplot_I_b,'String',iplot_I_b_s);
full_pipe_height_s=num2str(full_pipe_height); set(handles.full_pipe_height,'String',full_pipe_height_s);
%iRF_ele_s=num2str(iRF_ele); set(handles.iRF_ele,'String',iRF_ele_s);

function energy_Callback(hObject, eventdata, handles)
% hObject    handle to energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of energy as text
%        str2double(get(hObject,'String')) returns contents of energy as a double


% --- Executes during object creation, after setting all properties.
function energy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function I_b_Callback(hObject, eventdata, handles)
% hObject    handle to I_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of I_b as text
%        str2double(get(hObject,'String')) returns contents of I_b as a double


% --- Executes during object creation, after setting all properties.
function I_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to I_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function C_factor_Callback(hObject, eventdata, handles)
% hObject    handle to C_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of C_factor as text
%        str2double(get(hObject,'String')) returns contents of C_factor as a double


% --- Executes during object creation, after setting all properties.
function C_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function emit_norm_x_Callback(hObject, eventdata, handles)
% hObject    handle to emit_norm_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of emit_norm_x as text
%        str2double(get(hObject,'String')) returns contents of emit_norm_x as a double


% --- Executes during object creation, after setting all properties.
function emit_norm_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to emit_norm_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function emit_norm_y_Callback(hObject, eventdata, handles)
% hObject    handle to emit_norm_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of emit_norm_y as text
%        str2double(get(hObject,'String')) returns contents of emit_norm_y as a double


% --- Executes during object creation, after setting all properties.
function emit_norm_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to emit_norm_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_delta_Callback(hObject, eventdata, handles)
% hObject    handle to sigma_delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma_delta as text
%        str2double(get(hObject,'String')) returns contents of sigma_delta as a double


% --- Executes during object creation, after setting all properties.
function sigma_delta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma_delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function betax0_Callback(hObject, eventdata, handles)
% hObject    handle to betax0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of betax0 as text
%        str2double(get(hObject,'String')) returns contents of betax0 as a double


% --- Executes during object creation, after setting all properties.
function betax0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to betax0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function betay0_Callback(hObject, eventdata, handles)
% hObject    handle to betay0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of betay0 as text
%        str2double(get(hObject,'String')) returns contents of betay0 as a double


% --- Executes during object creation, after setting all properties.
function betay0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to betay0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alphax0_Callback(hObject, eventdata, handles)
% hObject    handle to alphax0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alphax0 as text
%        str2double(get(hObject,'String')) returns contents of alphax0 as a double


% --- Executes during object creation, after setting all properties.
function alphax0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphax0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alphay0_Callback(hObject, eventdata, handles)
% hObject    handle to alphay0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alphay0 as text
%        str2double(get(hObject,'String')) returns contents of alphay0 as a double


% --- Executes during object creation, after setting all properties.
function alphay0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphay0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function chirp_Callback(hObject, eventdata, handles)
% hObject    handle to chirp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chirp as text
%        str2double(get(hObject,'String')) returns contents of chirp as a double


% --- Executes during object creation, after setting all properties.
function chirp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chirp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function start_pos_Callback(hObject, eventdata, handles)
% hObject    handle to start_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_pos as text
%        str2double(get(hObject,'String')) returns contents of start_pos as a double


% --- Executes during object creation, after setting all properties.
function start_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function end_pos_Callback(hObject, eventdata, handles)
% hObject    handle to end_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_pos as text
%        str2double(get(hObject,'String')) returns contents of end_pos as a double


% --- Executes during object creation, after setting all properties.
function end_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rhox_Callback(hObject, eventdata, handles)
% hObject    handle to rhox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rhox as text
%        str2double(get(hObject,'String')) returns contents of rhox as a double


% --- Executes during object creation, after setting all properties.
function rhox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rhox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rhoy_Callback(hObject, eventdata, handles)
% hObject    handle to rhoy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rhoy as text
%        str2double(get(hObject,'String')) returns contents of rhoy as a double


% --- Executes during object creation, after setting all properties.
function rhoy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rhoy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ilattice.
function ilattice_Callback(hObject, eventdata, handles)
% hObject    handle to ilattice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ilattice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ilattice


% --- Executes during object creation, after setting all properties.
function ilattice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ilattice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over ilattice.
function ilattice_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ilattice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function lambda_start01_Callback(hObject, eventdata, handles)
% hObject    handle to lambda_start01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_start01 as text
%        str2double(get(hObject,'String')) returns contents of lambda_start01 as a double


% --- Executes during object creation, after setting all properties.
function lambda_start01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda_start01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lambda_end01_Callback(hObject, eventdata, handles)
% hObject    handle to lambda_end01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_end01 as text
%        str2double(get(hObject,'String')) returns contents of lambda_end01 as a double


% --- Executes during object creation, after setting all properties.
function lambda_end01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda_end01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scan_num01_Callback(hObject, eventdata, handles)
% hObject    handle to scan_num01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scan_num01 as text
%        str2double(get(hObject,'String')) returns contents of scan_num01 as a double


% --- Executes during object creation, after setting all properties.
function scan_num01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scan_num01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lambda_start02_Callback(hObject, eventdata, handles)
% hObject    handle to lambda_start02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_start02 as text
%        str2double(get(hObject,'String')) returns contents of lambda_start02 as a double


% --- Executes during object creation, after setting all properties.
function lambda_start02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda_start02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lambda_end02_Callback(hObject, eventdata, handles)
% hObject    handle to lambda_end02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_end02 as text
%        str2double(get(hObject,'String')) returns contents of lambda_end02 as a double


% --- Executes during object creation, after setting all properties.
function lambda_end02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda_end02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scan_num02_Callback(hObject, eventdata, handles)
% hObject    handle to scan_num02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scan_num02 as text
%        str2double(get(hObject,'String')) returns contents of scan_num02 as a double


% --- Executes during object creation, after setting all properties.
function scan_num02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scan_num02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mesh_num_Callback(hObject, eventdata, handles)
% hObject    handle to mesh_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mesh_num as text
%        str2double(get(hObject,'String')) returns contents of mesh_num as a double


% --- Executes during object creation, after setting all properties.
function mesh_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mesh_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iplot_lattice_Callback(hObject, eventdata, handles)
% hObject    handle to iplot_lattice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iplot_lattice as text
%        str2double(get(hObject,'String')) returns contents of iplot_lattice as a double


% --- Executes during object creation, after setting all properties.
function iplot_lattice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iplot_lattice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iplot_gain_func_Callback(hObject, eventdata, handles)
% hObject    handle to iplot_gain_func (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iplot_gain_func as text
%        str2double(get(hObject,'String')) returns contents of iplot_gain_func as a double


% --- Executes during object creation, after setting all properties.
function iplot_gain_func_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iplot_gain_func (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iplot_gain_spec_Callback(hObject, eventdata, handles)
% hObject    handle to iplot_gain_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iplot_gain_spec as text
%        str2double(get(hObject,'String')) returns contents of iplot_gain_spec as a double


% --- Executes during object creation, after setting all properties.
function iplot_gain_spec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iplot_gain_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function isave_Callback(hObject, eventdata, handles)
% hObject    handle to isave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of isave as text
%        str2double(get(hObject,'String')) returns contents of isave as a double


% --- Executes during object creation, after setting all properties.
function isave_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ilast_Callback(hObject, eventdata, handles)
% hObject    handle to ilast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ilast as text
%        str2double(get(hObject,'String')) returns contents of ilast as a double


% --- Executes during object creation, after setting all properties.
function ilast_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ilast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iiterative_Callback(hObject, eventdata, handles)
% hObject    handle to iiterative (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iiterative as text
%        str2double(get(hObject,'String')) returns contents of iiterative as a double


% --- Executes during object creation, after setting all properties.
function iiterative_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iiterative (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iCSR_ss_Callback(hObject, eventdata, handles)
% hObject    handle to iCSR_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iCSR_ss as text
%        str2double(get(hObject,'String')) returns contents of iCSR_ss as a double


% --- Executes during object creation, after setting all properties.
function iCSR_ss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iCSR_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iCSR_tr_Callback(hObject, eventdata, handles)
% hObject    handle to iCSR_tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iCSR_tr as text
%        str2double(get(hObject,'String')) returns contents of iCSR_tr as a double


% --- Executes during object creation, after setting all properties.
function iCSR_tr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iCSR_tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iCSR_drift_Callback(hObject, eventdata, handles)
% hObject    handle to iCSR_drift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iCSR_drift as text
%        str2double(get(hObject,'String')) returns contents of iCSR_drift as a double


% --- Executes during object creation, after setting all properties.
function iCSR_drift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iCSR_drift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iLSC_Callback(hObject, eventdata, handles)
% hObject    handle to iLSC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iLSC as text
%        str2double(get(hObject,'String')) returns contents of iLSC as a double


% --- Executes during object creation, after setting all properties.
function iLSC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iLSC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ilinac_Callback(hObject, eventdata, handles)
% hObject    handle to ilinac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ilinac as text
%        str2double(get(hObject,'String')) returns contents of ilinac as a double


% --- Executes during object creation, after setting all properties.
function ilinac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ilinac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iGaussian_Callback(hObject, eventdata, handles)
% hObject    handle to iGaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iGaussian as text
%        str2double(get(hObject,'String')) returns contents of iGaussian as a double


% --- Executes during object creation, after setting all properties.
function iGaussian_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iGaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iEnergy_mod_calc_Callback(hObject, eventdata, handles)
% hObject    handle to iEnergy_mod_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iEnergy_mod_calc as text
%        str2double(get(hObject,'String')) returns contents of iEnergy_mod_calc as a double


% --- Executes during object creation, after setting all properties.
function iEnergy_mod_calc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iEnergy_mod_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iplot_energy_mod_Callback(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iplot_energy_mod as text
%        str2double(get(hObject,'String')) returns contents of iplot_energy_mod as a double


% --- Executes during object creation, after setting all properties.
function iplot_energy_mod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%----------------
function idm_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iplot_energy_mod as text
%        str2double(get(hObject,'String')) returns contents of iplot_energy_mod as a double


% --- Executes during object creation, after setting all properties.
function idm_analysis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function surfplot_gain_spec_func_Callback(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iplot_energy_mod as text
%        str2double(get(hObject,'String')) returns contents of iplot_energy_mod as a double


% --- Executes during object creation, after setting all properties.
function surfplot_gain_spec_func_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function issCSRpp_Callback(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iplot_energy_mod as text
%        str2double(get(hObject,'String')) returns contents of iplot_energy_mod as a double


% --- Executes during object creation, after setting all properties.
function issCSRpp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ssCSR_model_Callback(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iplot_energy_mod as text
%        str2double(get(hObject,'String')) returns contents of iplot_energy_mod as a double


% --- Executes during object creation, after setting all properties.
function ssCSR_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LSC_model_Callback(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iplot_energy_mod as text
%        str2double(get(hObject,'String')) returns contents of iplot_energy_mod as a double


% --- Executes during object creation, after setting all properties.
function LSC_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function iquilt_plot_Callback(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iplot_energy_mod as text
%        str2double(get(hObject,'String')) returns contents of iplot_energy_mod as a double


% --- Executes during object creation, after setting all properties.
function iquilt_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function iplot_I_b_Callback(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iplot_energy_mod as text
%        str2double(get(hObject,'String')) returns contents of iplot_energy_mod as a double


% --- Executes during object creation, after setting all properties.
function iplot_I_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function full_pipe_height_Callback(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iplot_energy_mod as text
%        str2double(get(hObject,'String')) returns contents of iplot_energy_mod as a double


% --- Executes during object creation, after setting all properties.
function full_pipe_height_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function iRF_ele_Callback(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iplot_energy_mod as text
%        str2double(get(hObject,'String')) returns contents of iplot_energy_mod as a double


% --- Executes during object creation, after setting all properties.
function iRF_ele_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iplot_energy_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
%axes(handles.axes1)
%matlabImage = imread('lambda_scan_illustration.png');
%image(matlabImage)
imshow('lambda_scan_illustration.png');
axis off
axis image



function iEnergy_gain_calc_Callback(hObject, eventdata, handles)
% hObject    handle to iEnergy_gain_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iEnergy_gain_calc as text
%        str2double(get(hObject,'String')) returns contents of iEnergy_gain_calc as a double



% --- Executes during object creation, after setting all properties.
function iEnergy_gain_calc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iEnergy_gain_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when uipanel2 is resized.
function uipanel2_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in iIBS.
function iIBS_Callback(hObject, eventdata, handles)
% hObject    handle to iIBS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of iIBS


% --- Executes on button press in DO_X.
function DO_X_Callback(hObject, eventdata, handles)
% hObject    handle to DO_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DO_X


% --- Executes on button press in DO_Y.
function DO_Y_Callback(hObject, eventdata, handles)
% hObject    handle to DO_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DO_Y


% --- Executes on button press in DO_Z.
function DO_Z_Callback(hObject, eventdata, handles)
% hObject    handle to DO_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DO_Z


% --- Executes on button press in iDz.
function iDz_Callback(hObject, eventdata, handles)
% hObject    handle to iDz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of iDz


% --- Executes on button press in iDzz.
function iDzz_Callback(hObject, eventdata, handles)
% hObject    handle to iDzz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of iDzz


% --- Executes on button press in iSES.
function iSES_Callback(hObject, eventdata, handles)
% hObject    handle to iSES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of iSES


% --- Executes on button press in itransLD.
function itransLD_Callback(hObject, eventdata, handles)
% hObject    handle to itransLD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of itransLD


% --- Executes on button press in plot_Zk_vs_s.
function plot_Zk_vs_s_Callback(hObject, eventdata, handles)
% hObject    handle to plot_Zk_vs_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_Zk_vs_s


% --- Executes on button press in use_orig_mat_output.
function use_orig_mat_output_Callback(hObject, eventdata, handles)
% hObject    handle to use_orig_mat_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_orig_mat_output



function CSRdrifmodel_Callback(hObject, eventdata, handles)
% hObject    handle to CSRdrifmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CSRdrifmodel as text
%        str2double(get(hObject,'String')) returns contents of CSRdrifmodel as a double


% --- Executes during object creation, after setting all properties.
function CSRdrifmodel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CSRdrifmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in iFELcalc.
function iFELcalc_Callback(hObject, eventdata, handles)
% hObject    handle to iFELcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of iFELcalc



function edit52_Callback(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit52 as text
%        str2double(get(hObject,'String')) returns contents of edit52 as a double


% --- Executes during object creation, after setting all properties.
function edit52_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Derbenev.
function Derbenev_Callback(hObject, eventdata, handles)
% hObject    handle to Derbenev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Derbenev


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8


% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox31.
function checkbox31_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox31


% --- Executes on selection change in gen_CSR_LSC_param.
function gen_CSR_LSC_param_Callback(hObject, eventdata, handles)
% hObject    handle to gen_CSR_LSC_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns gen_CSR_LSC_param contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gen_CSR_LSC_param


% --- Executes during object creation, after setting all properties.
function gen_CSR_LSC_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gen_CSR_LSC_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox32.
function checkbox32_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox32
