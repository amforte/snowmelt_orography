% Requires 
% API key from OpenTopography that is input as api_key='string';
% TopoToolbox and Topographic-Analysis-Kit (and the required toolboxes for those)

GC_DEM=readopentopo('north',45,'south',39.5,'west',38,'east',51,'demtype','SRTMGL3','apikey',api_key);
ALPS_DEM=readopentopo('north',50,'south',43,'west',5,'east',16,'demtype','SRTMGL3','apikey',api_key);
BC_DEM=readopentopo('north',54,'south',48,'west',-131,'east',-120,'demtype','SRTMGL3','apikey',api_key);

GC_DEM=reproject2utm(GC_DEM,90);
ALPS_DEM=reproject2utm(ALPS_DEM,90);
BC_DEM=reproject2utm(BC_DEM,90);

[DEM,FD,A,S]=MakeStreams(GC_DEM,1e6,'file_name','GC','no_data_exp','DEM<=250','mex',true);
[DEM,FD,A,S]=MakeStreams(ALPS_DEM,1e6,'file_name','ALPS','no_data_exp','DEM<=250','mex',true);
[DEM,FD,A,S]=MakeStreams(BC_DEM,1e6,'file_name','BC','no_data_exp','DEM<=250','mex',true);

save('GC.mat','GC_DEM','-append');
save('ALPS.mat','ALPS_DEM','-append');
save('BC.mat','BC_DEM','-append');

load('GC.mat','DEM','FD','A','S');
ProcessRiverBasins(DEM,FD,A,S,300,'gc_basins','calc_relief',true,'min_basin_size',50);
SubDivideBigBasins('gc_basins',250,'filtered_trunk','SBFiles_Dir','sub250','min_basin_size',25);

load('ALPS.mat','DEM','FD','A','S');
ProcessRiverBasins(DEM,FD,A,S,300,'alps_basins','calc_relief',true,'min_basin_size',50);
SubDivideBigBasins('alps_basins',250,'filtered_trunk','SBFiles_Dir','sub250','min_basin_size',25);

load('BC.mat','DEM','FD','A','S');
ProcessRiverBasins(DEM,FD,A,S,300,'bc_basins','calc_relief',true,'min_basin_size',50);
SubDivideBigBasins('bc_basins',250,'filtered_trunk','SBFiles_Dir','sub250','min_basin_size',25);

[T_GC]=CompileBasinStats('gc_basins','location_of_subbasins','sub250','include','subdivided');
[T_ALPS]=CompileBasinStats('alps_basins','location_of_subbasins','sub250','include','subdivided');
[T_BC]=CompileBasinStats('bc_basins','location_of_subbasins','sub250','include','subdivided');

writetable(T_GC,'gc_ksn_rlf.csv');
writetable(T_ALPS,'alps_ksn_rlf.csv');
writetable(T_BC,'bc_ksn_rlf.csv');