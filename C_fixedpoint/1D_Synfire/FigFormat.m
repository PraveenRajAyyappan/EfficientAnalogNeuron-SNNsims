function FigFormat(type)
%% FigFormat.m::Function to format figures into standardized sizes
% Written Nov 2021 by: Swagat Bhattacharyya
% Notes: type==1 => Portrait mode, type==2 => Landscape mode

%% Get figure handle
Handle = gcf;

%% Set default line widths
box on;
set(Handle,'defaultlinelinewidth',2)
set(Handle,'defaultaxeslinewidth',1)
set(findobj(gca,'Type','line'),'linewidth',2)
set(findobj(gca,'Type','axes'),'linewidth',1)

%% Resize figure
if type==1
    set(Handle, 'Units','centimeters', 'Position',[8 0 16 16])
else
    set(Handle, 'Units','centimeters', 'Position',[8 6 16.18 10])
end

%% Set color scheme
set(Handle,'color','w')

%% Set font sizes
set(gca, 'Fontsize',12)
title(get(get(gca,'Title'),'String'),'Fontsize',12)

% Set properties to make any export look clean
set(gcf, 'PaperPositionMode', 'auto');
end