function varargout = GUI_FEL_calc(varargin)
% GUI_FEL_CALC MATLAB code for GUI_FEL_calc.fig
%      GUI_FEL_CALC, by itself, creates a new GUI_FEL_CALC or raises the existing
%      singleton*.
%
%      H = GUI_FEL_CALC returns the handle to a new GUI_FEL_CALC or the handle to
%      the existing singleton*.
%
%      GUI_FEL_CALC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FEL_CALC.M with the given input arguments.
%
%      GUI_FEL_CALC('Property','Value',...) creates a new GUI_FEL_CALC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_FEL_calc_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_FEL_calc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_FEL_calc

% Last Modified by GUIDE v2.5 23-Aug-2022 13:15:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_FEL_calc_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_FEL_calc_OutputFcn, ...
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


% --- Executes just before GUI_FEL_calc is made visible.
function GUI_FEL_calc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_FEL_calc (see VARARGIN)

% Choose default command line output for GUI_FEL_calc
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_FEL_calc wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_FEL_calc_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

global MX und_lambda und_L lambda_r P_L bk_ini Pk_ini
global I_b_fin sigma_delta_fin sig_x emit_nx_fin und_beta_ave

I_b_fin_s=num2str(I_b_fin); set(handles.I_b_fin,'String',I_b_fin_s);
sigma_delta_fin_s=num2str(sigma_delta_fin); set(handles.sigma_delta_fin,'String',sigma_delta_fin_s);
sig_x_s=num2str(sig_x); set(handles.sig_x,'String',sig_x_s);
emit_nx_fin_s=num2str(emit_nx_fin); set(handles.emit_nx_fin,'String',emit_nx_fin_s);
bk_ini_s=num2str(bk_ini); set(handles.bk_ini,'String',bk_ini_s);
Pk_ini_s=num2str(Pk_ini); set(handles.Pk_ini,'String',Pk_ini_s);
und_beta_ave_s=num2str(und_beta_ave); set(handles.und_beta_ave,'String',und_beta_ave_s);
und_lambda_s=num2str(und_lambda); set(handles.und_lambda,'String',und_lambda_s);
und_L_s=num2str(und_L); set(handles.und_L,'String',und_L_s);
lambda_r_s=num2str(lambda_r); set(handles.lambda_r,'String',lambda_r_s);
P_L_s=num2str(P_L); set(handles.P_L,'String',P_L_s);
MX=handles.MX.Value;



function und_beta_ave_Callback(hObject, eventdata, handles)
% hObject    handle to und_beta_ave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of und_beta_ave as text
%        str2double(get(hObject,'String')) returns contents of und_beta_ave as a double


% --- Executes during object creation, after setting all properties.
function und_beta_ave_CreateFcn(hObject, eventdata, handles)
% hObject    handle to und_beta_ave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function und_lambda_Callback(hObject, eventdata, handles)
% hObject    handle to und_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of und_lambda as text
%        str2double(get(hObject,'String')) returns contents of und_lambda as a double


% --- Executes during object creation, after setting all properties.
function und_lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to und_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function und_L_Callback(hObject, eventdata, handles)
% hObject    handle to und_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of und_L as text
%        str2double(get(hObject,'String')) returns contents of und_L as a double


% --- Executes during object creation, after setting all properties.
function und_L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to und_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lambda_r_Callback(hObject, eventdata, handles)
% hObject    handle to lambda_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_r as text
%        str2double(get(hObject,'String')) returns contents of lambda_r as a double


% --- Executes during object creation, after setting all properties.
function lambda_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function I_b_fin_Callback(hObject, eventdata, handles)
% hObject    handle to I_b_fin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of I_b_fin as text
%        str2double(get(hObject,'String')) returns contents of I_b_fin as a double


% --- Executes during object creation, after setting all properties.
function I_b_fin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to I_b_fin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_delta_fin_Callback(hObject, eventdata, handles)
% hObject    handle to sigma_delta_fin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma_delta_fin as text
%        str2double(get(hObject,'String')) returns contents of sigma_delta_fin as a double


% --- Executes during object creation, after setting all properties.
function sigma_delta_fin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma_delta_fin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function emit_nx_fin_Callback(hObject, eventdata, handles)
% hObject    handle to emit_nx_fin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of emit_nx_fin as text
%        str2double(get(hObject,'String')) returns contents of emit_nx_fin as a double


% --- Executes during object creation, after setting all properties.
function emit_nx_fin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to emit_nx_fin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sig_x_Callback(hObject, eventdata, handles)
% hObject    handle to sig_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sig_x as text
%        str2double(get(hObject,'String')) returns contents of sig_x as a double


% --- Executes during object creation, after setting all properties.
function sig_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sig_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MX.
function MX_Callback(hObject, eventdata, handles)
% hObject    handle to MX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MX



function P_L_Callback(hObject, eventdata, handles)
% hObject    handle to P_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of P_L as text
%        str2double(get(hObject,'String')) returns contents of P_L as a double


% --- Executes during object creation, after setting all properties.
function P_L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to P_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bk_ini_Callback(hObject, eventdata, handles)
% hObject    handle to bk_ini (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bk_ini as text
%        str2double(get(hObject,'String')) returns contents of bk_ini as a double


% --- Executes during object creation, after setting all properties.
function bk_ini_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bk_ini (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pk_ini_Callback(hObject, eventdata, handles)
% hObject    handle to Pk_ini (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pk_ini as text
%        str2double(get(hObject,'String')) returns contents of Pk_ini as a double


% --- Executes during object creation, after setting all properties.
function Pk_ini_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pk_ini (see GCBO)
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
