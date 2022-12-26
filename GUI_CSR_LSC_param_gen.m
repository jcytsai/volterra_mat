function varargout = GUI_CSR_LSC_param_gen(varargin)
% GUI_CSR_LSC_PARAM_GEN MATLAB code for GUI_CSR_LSC_param_gen.fig
%      GUI_CSR_LSC_PARAM_GEN, by itself, creates a new GUI_CSR_LSC_PARAM_GEN or raises the existing
%      singleton*.
%
%      H = GUI_CSR_LSC_PARAM_GEN returns the handle to a new GUI_CSR_LSC_PARAM_GEN or the handle to
%      the existing singleton*.
%
%      GUI_CSR_LSC_PARAM_GEN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CSR_LSC_PARAM_GEN.M with the given input arguments.
%
%      GUI_CSR_LSC_PARAM_GEN('Property','Value',...) creates a new GUI_CSR_LSC_PARAM_GEN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_CSR_LSC_param_gen_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_CSR_LSC_param_gen_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_CSR_LSC_param_gen

% Last Modified by GUIDE v2.5 26-Aug-2022 12:44:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_CSR_LSC_param_gen_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_CSR_LSC_param_gen_OutputFcn, ...
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


% --- Executes just before GUI_CSR_LSC_param_gen is made visible.
function GUI_CSR_LSC_param_gen_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_CSR_LSC_param_gen (see VARARGIN)

% Choose default command line output for GUI_CSR_LSC_param_gen
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_CSR_LSC_param_gen wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_CSR_LSC_param_gen_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%global Gf_max lambda_opt lambda_min num_bin_per_lambda num_ptc_per_bin sig_s_ini 
global num_kick num_bin num_ptc cutoff_01 cutoff_02 ini_dmod_amp ini_emod_amp

num_ptc_s=num2str(num_ptc); set(handles.num_ptc,'String',num_ptc_s);
num_kick_s=num2str(num_kick); set(handles.num_kick,'String',num_kick_s);
num_bin_s=num2str(num_bin); set(handles.num_bin,'String',num_bin_s);
cutoff_01_s=num2str(cutoff_01); set(handles.cutoff_01,'String',cutoff_01_s);
cutoff_02_s=num2str(cutoff_02); set(handles.cutoff_02,'String',cutoff_02_s);
ini_dmod_amp_s=num2str(ini_dmod_amp); set(handles.ini_dmod_amp,'String',ini_dmod_amp_s);


function num_ptc_Callback(hObject, eventdata, handles)
% hObject    handle to num_ptc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_ptc as text
%        str2double(get(hObject,'String')) returns contents of num_ptc as a double


% --- Executes during object creation, after setting all properties.
function num_ptc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_ptc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_kick_Callback(hObject, eventdata, handles)
% hObject    handle to num_kick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_kick as text
%        str2double(get(hObject,'String')) returns contents of num_kick as a double


% --- Executes during object creation, after setting all properties.
function num_kick_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_kick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_bin_Callback(hObject, eventdata, handles)
% hObject    handle to num_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_bin as text
%        str2double(get(hObject,'String')) returns contents of num_bin as a double


% --- Executes during object creation, after setting all properties.
function num_bin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cutoff_01_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff_01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cutoff_01 as text
%        str2double(get(hObject,'String')) returns contents of cutoff_01 as a double


% --- Executes during object creation, after setting all properties.
function cutoff_01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoff_01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cutoff_02_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cutoff_02 as text
%        str2double(get(hObject,'String')) returns contents of cutoff_02 as a double


% --- Executes during object creation, after setting all properties.
function cutoff_02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoff_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ini_dmod_amp_Callback(hObject, eventdata, handles)
% hObject    handle to ini_dmod_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ini_dmod_amp as text
%        str2double(get(hObject,'String')) returns contents of ini_dmod_amp as a double


% --- Executes during object creation, after setting all properties.
function ini_dmod_amp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ini_dmod_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to num_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_bin as text
%        str2double(get(hObject,'String')) returns contents of num_bin as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff_01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cutoff_01 as text
%        str2double(get(hObject,'String')) returns contents of cutoff_01 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoff_01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cutoff_02 as text
%        str2double(get(hObject,'String')) returns contents of cutoff_02 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoff_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
