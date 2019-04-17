 function varargout = MPM_FD_FLIM_SNR(varargin)
% MPM_FD_FLIM_SNR MATLAB code for MPM_FD_FLIM_SNR.fig
%
%   Author: Yide Zhang
%   Email: yzhang34@nd.edu
%   Date: April 16, 2019
%   Copyright: University of Notre Dame, 2019

% Last Modified by GUIDE v2.5 16-Apr-2019 21:13:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MPM_FD_FLIM_SNR_OpeningFcn, ...
                   'gui_OutputFcn',  @MPM_FD_FLIM_SNR_OutputFcn, ...
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


% --- Executes just before MPM_FD_FLIM_SNR is made visible.
function MPM_FD_FLIM_SNR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MPM_FD_FLIM_SNR (see VARARGIN)

% Choose default command line output for MPM_FD_FLIM_SNR
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MPM_FD_FLIM_SNR wait for user response (see UIRESUME)
% uiwait(handles.figure1);
clc

addpath('./functions')

set(hObject,'Toolbar','figure'); % let the toolbar be operable
plot_MPM_FLIM_Phase(handles, false)


% --- Outputs from this function are returned to the command line.
function varargout = MPM_FD_FLIM_SNR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Popup_Exc.
function Popup_Exc_Callback(hObject, eventdata, handles)
% hObject    handle to Popup_Exc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Popup_Exc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Popup_Exc
plot_MPM_FLIM_Phase(handles, false)


% --- Executes during object creation, after setting all properties.
function Popup_Exc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Popup_Exc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditText_f_mod_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_f_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_f_mod as text
%        str2double(get(hObject,'String')) returns contents of EditText_f_mod as a double
f_mod = str2double(get(hObject,'String'));
T_mod = 1/f_mod;
set(handles.EditText_T_mod, 'String', num2str(T_mod,'%10.2e'));
plot_MPM_FLIM_Phase(handles, false)



% --- Executes during object creation, after setting all properties.
function EditText_f_mod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_f_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditText_f_sample_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_f_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_f_sample as text
%        str2double(get(hObject,'String')) returns contents of EditText_f_sample as a double
f_sample = str2double(get(hObject,'String'));
T_sample = 1/f_sample;
set(handles.EditText_T_sample, 'String', num2str(T_sample,'%10.2e'));
plot_MPM_FLIM_Phase(handles, false)


% --- Executes during object creation, after setting all properties.
function EditText_f_sample_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_f_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditText_N_period_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_N_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_N_period as text
%        str2double(get(hObject,'String')) returns contents of EditText_N_period as a double
plot_MPM_FLIM_Phase(handles, false)


% --- Executes during object creation, after setting all properties.
function EditText_N_period_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_N_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditText_T_mod_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_T_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_T_mod as text
%        str2double(get(hObject,'String')) returns contents of EditText_T_mod as a double
T_mod = str2double(get(hObject,'String'));
f_mod = 1/T_mod;
set(handles.EditText_f_mod, 'String', num2str(f_mod,'%10.2e'));
plot_MPM_FLIM_Phase(handles, false)


% --- Executes during object creation, after setting all properties.
function EditText_T_mod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_T_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function EditText_T_sample_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_T_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_T_sample as text
%        str2double(get(hObject,'String')) returns contents of EditText_T_sample as a double
T_sample = str2double(get(hObject,'String'));
f_sample = 1/T_sample;
set(handles.EditText_f_sample, 'String', num2str(f_sample,'%10.2e'));
plot_MPM_FLIM_Phase(handles, false)


% --- Executes during object creation, after setting all properties.
function EditText_T_sample_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_T_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditText_N_frame_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_N_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_N_frame as text
%        str2double(get(hObject,'String')) returns contents of EditText_N_frame as a double


% --- Executes during object creation, after setting all properties.
function EditText_N_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_N_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditText_phi_shift_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_phi_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_phi_shift as text
%        str2double(get(hObject,'String')) returns contents of EditText_phi_shift as a double
plot_MPM_FLIM_Phase(handles, false)


% --- Executes during object creation, after setting all properties.
function EditText_phi_shift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_phi_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Button_Calc.
function Button_Calc_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_MPM_FLIM_Phase(handles, true)



function EditText_m_DC_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_m_DC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_m_DC as text
%        str2double(get(hObject,'String')) returns contents of EditText_m_DC as a double
m_DC = str2double(get(hObject,'String'));
m_AC = str2double(get(handles.EditText_m_AC,'String'));
m = m_AC/m_DC;
set(handles.EditText_m, 'String', num2str(m,'%10.2f'));
plot_MPM_FLIM_Phase(handles, false)


% --- Executes during object creation, after setting all properties.
function EditText_m_DC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_m_DC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function EditText_m_AC_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_m_AC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_m_AC as text
%        str2double(get(hObject,'String')) returns contents of EditText_m_AC as a double
m_AC = str2double(get(hObject,'String'));
m_DC = str2double(get(handles.EditText_m_DC,'String'));
m = m_AC/m_DC;
set(handles.EditText_m, 'String', num2str(m,'%10.2f'));
plot_MPM_FLIM_Phase(handles, false)


% --- Executes during object creation, after setting all properties.
function EditText_m_AC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_m_AC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditText_m_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_m as text
%        str2double(get(hObject,'String')) returns contents of EditText_m as a double
m = str2double(get(hObject,'String'));
m_DC = str2double(get(handles.EditText_m_DC,'String'));
m_AC = m*m_DC;
set(handles.EditText_m_AC, 'String', num2str(m_AC,'%10.2f'));
plot_MPM_FLIM_Phase(handles, false)


% --- Executes during object creation, after setting all properties.
function EditText_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditText_a_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_a as text
%        str2double(get(hObject,'String')) returns contents of EditText_a as a double
plot_MPM_FLIM_Phase(handles, false)


% --- Executes during object creation, after setting all properties.
function EditText_a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditText_real_tau_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_real_tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_real_tau as text
%        str2double(get(hObject,'String')) returns contents of EditText_real_tau as a double
plot_MPM_FLIM_Phase(handles, false)


% --- Executes during object creation, after setting all properties.
function EditText_real_tau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_real_tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditText_eff_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_eff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_eff as text
%        str2double(get(hObject,'String')) returns contents of EditText_eff as a double
plot_MPM_FLIM_Phase(handles, false)


% --- Executes during object creation, after setting all properties.
function EditText_eff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_eff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditText_psn_noise_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_psn_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_psn_noise as text
%        str2double(get(hObject,'String')) returns contents of EditText_psn_noise as a double
plot_MPM_FLIM_Phase(handles, false)


% --- Executes during object creation, after setting all properties.
function EditText_psn_noise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_psn_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Popup_Frame.
function Popup_Frame_Callback(hObject, eventdata, handles)
% hObject    handle to Popup_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Popup_Frame contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Popup_Frame

% --- Executes during object creation, after setting all properties.
function Popup_Frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Popup_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Check_Log.
function Check_Log_Callback(hObject, eventdata, handles)
% hObject    handle to Check_Log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Check_Log



function EditText_N_int_Callback(hObject, eventdata, handles)
% hObject    handle to EditText_N_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditText_N_int as text
%        str2double(get(hObject,'String')) returns contents of EditText_N_int as a double


% --- Executes during object creation, after setting all properties.
function EditText_N_int_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditText_N_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Check_Draw.
function Check_Draw_Callback(hObject, eventdata, handles)
% hObject    handle to Check_Draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Check_Draw
