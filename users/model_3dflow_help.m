function varargout = model_3dflow_help(varargin)
% MODEL_3DFLOW_HELP MATLAB code for model_3dflow_help.fig
%      MODEL_3DFLOW_HELP, by itself, creates a new MODEL_3DFLOW_HELP or raises the existing
%      singleton*.
%
%      H = MODEL_3DFLOW_HELP returns the handle to a new MODEL_3DFLOW_HELP or the handle to
%      the existing singleton*.
%
%      MODEL_3DFLOW_HELP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODEL_3DFLOW_HELP.M with the given input arguments.
%
%      MODEL_3DFLOW_HELP('Property','Value',...) creates a new MODEL_3DFLOW_HELP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before model_3dflow_help_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to model_3dflow_help_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help model_3dflow_help

% Last Modified by GUIDE v2.5 17-Dec-2021 16:41:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @model_3dflow_help_OpeningFcn, ...
                   'gui_OutputFcn',  @model_3dflow_help_OutputFcn, ...
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


% --- Executes just before model_3dflow_help is made visible.
function model_3dflow_help_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to model_3dflow_help (see VARARGIN)

% Choose default command line output for model_3dflow_help
handles.output = hObject;

if varargin{1} == 1
    set(handles.help3d_topo,'Visible','on');
    set(handles.help3d_init,'Visible','off');
    set(handles.help3d_params,'Visible','off');
elseif varargin{1} == 2
    set(handles.help3d_topo,'Visible','off');
    set(handles.help3d_init,'Visible','on');
    set(handles.help3d_params,'Visible','off');
else
    set(handles.help3d_topo,'Visible','off');
    set(handles.help3d_init,'Visible','off');
    set(handles.help3d_params,'Visible','on');
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes model_3dflow_help wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = model_3dflow_help_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
