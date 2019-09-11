function var=ncreadall(filename,varargin)
% This is a function to fully read a NetCDF file using matlab version
% later than 2011
%% Input
% filename: name of the NetCDF file.
% varargin: consist of 2 parts. The first part is [{start},{count},{stride}]
%           the second part is the serial number that you want to specify [{start},{count},{stride}]. 
%% Output
% var: struct array containing all the variables that have been read,
% having been applied add_offset, scale_factor and missing_data
%% Examples
%% Copyright
% By Hong Shen, Dec. 01, 2013, catDoraemon1991@gmail.com
% Open source code, please contact the author when any change is made
finfo=ncinfo(filename);
if nargin==1
  for i=1:numel(finfo.Variables)
      varname=finfo.Variables(i).Name;
      if isempty(str2num(varname(1)))
      var.(varname)=ncread(filename,varname);
      end
  end
else
    if varargin{end}>numel(finfo.Variables)||varargin{end}<1
        error(['Please input a integer between 1 and ' num2str(numel(finfo.Variables)) ' for the variable you want to assign start, count or stride']);
    end
    for i=[1:varargin{end}-1,varargin{end}+1:numel(finfo.Variables)]
        varname=finfo.Variables(i).Name;
        var.(varname)=ncread(filename,varname);
    end
    varname=finfo.Variables(varargin{end}).Name;
    var.(varname)=ncread(filename,varname,varargin{1:end-1});
end
end

