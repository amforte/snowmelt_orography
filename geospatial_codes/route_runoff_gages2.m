function route_runoff_gages2(strt)
	% Written by Adam M. Forte
	% aforte8@lsu.edu
	%
	% This function downloads DEMs from OpenTopography that cover individual drainage basins
	% within a select portion of the HCDN-2009 subset of GagesII and then uses daily WaterGAP3
	% rasters to route runoff to the outlet of individual watersheds and builds an expected
	% daily discharge and runoff record for each watershed. 	
	%
	% Requires: 
	% An API key from OpenTopography that is input as api_key='string';
	% gagesII data, available here:
	% https://cmerwebmap.cr.usgs.gov/catalog/item/5788f619e4b0d27deb389055
	% Water Resources Reanalysis v2 - WaterGAP3 Data, available here:
	% http://www.earth2observe.eu/

	% Load shapefile
	shp=shaperead(['boundaries_shapefiles_by_aggeco' filesep 'bas_ref_all_wgs84.shp']);
	% Load in list of filtered HCDN-2009 basins
	hcdn=readtable('filtered_gages2_staid.csv');
	% Build lid of staid
	for ii=1:numel(shp)
		ref_staid(ii,1)=str2num(shp(ii,1).GAGE_ID);
	end
	% Find index of hcdn basins in shape
	idx=ismember(ref_staid,hcdn.STAID);
	shp=shp(idx,1);

	% Define file locations for nc files
	Qs1=['e2o_univk_wrr2' filesep 'e2o_univk_wrr2_glob15_day_Qs_1980-1989.nc'];
	Qs2=['e2o_univk_wrr2' filesep 'e2o_univk_wrr2_glob15_day_Qs_1990-1999.nc'];
	Qsb1=['e2o_univk_wrr2' filesep 'e2o_univk_wrr2_glob15_day_Qsb_1980-1989.nc'];
	Qsb2=['e2o_univk_wrr2' filesep 'e2o_univk_wrr2_glob15_day_Qsb_1990-1999.nc'];
	Qsm1=['e2o_univk_wrr2' filesep 'e2o_univk_wrr2_glob15_day_Qsm_1980-1989.nc'];
	Qsm2=['e2o_univk_wrr2' filesep 'e2o_univk_wrr2_glob15_day_Qsm_1990-1999.nc'];

	lat=ncread(Qs1,'lat');
	lon=ncread(Qs1,'lon');

	num_files=numel(shp);

	for ii=strt:num_files
		disp(['Working on ' shp(ii,1).GAGE_ID])

		% Extract bounding box
		bb=shp(ii,1).BoundingBox;
		lft=bb(1,1); rgt=bb(2,1); top=bb(2,2); bot=bb(1,2);
		shpc=shp(ii,1);

		DEMLL=readopentopo('north',top,'south',bot,'west',lft,'east',rgt,'demtype','SRTMGL3','apikey',api_key,'verbose',false);
		C=polygon2GRIDobj(DEMLL,shpc);
		DEMc=clip(DEMLL,C);
		DEM=reproject2utm(DEMc,90);
		FD=FLOWobj(DEM,'preprocess','carve','mex',true);
		DA=flowacc(FD)*(DEM.cellsize^2);
		S=STREAMobj(FD,'minarea',1e6,'unit','mapunits');
		oix=streampoi(klargestconncomps(S,1),'outlets','ix');

		% Find indices within nc
		[~,lft_ix]=min(abs(lon-lft));
		[~,rgt_ix]=min(abs(lon-rgt));
		[~,top_ix]=min(abs(lat-top));
		[~,bot_ix]=min(abs(lat-bot));

		% GRIDobj apparently can't handle a single column/row grid. Inflate by one if either
		% index is equal to it's adjacent side (the resampling will trim it to the proper dimensions)
		if lft_ix==rgt_ix
			rgt_ix=rgt_ix+1;
		end

		if top_ix==bot_ix
			top_ix=top_ix+1;
		end

		startL=[lft_ix bot_ix 1];
		countL=[(rgt_ix-lft_ix)+1 (top_ix-bot_ix)+1 inf];

		% Grab Time Series and convert to m/day
		Qs1m=(ncread(Qs1,'Qs',startL,countL)./1000)*(60*60*24)*-1;
		Qs2m=(ncread(Qs2,'Qs',startL,countL)./1000)*(60*60*24)*-1;
		Qsm=cat(3,Qs1m,Qs2m);
		disp('	Qs complete')
		Qsb1m=(ncread(Qsb1,'Qsb',startL,countL)./1000)*(60*60*24)*-1;
		Qsb2m=(ncread(Qsb2,'Qsb',startL,countL)./1000)*(60*60*24)*-1;
		Qsbm=cat(3,Qsb1m,Qsb2m);
		disp('	Qsb complete')
		Qsm1m=(ncread(Qsm1,'Qsm',startL,countL)./1000)*(60*60*24);
		Qsm2m=(ncread(Qsm2,'Qsm',startL,countL)./1000)*(60*60*24);
		Qsmm=cat(3,Qsm1m,Qsm2m);
		disp('	Qsm complete')

		num_days=size(Qsm,3);

		runoff=zeros(num_days,1);

		for jj=1:num_days
			% Grab day and transpose
			if mod(jj,100)==0
				disp(['	' num2str(jj) ' of ' num2str(num_days)])
			end
			Qs0=Qsm(:,:,jj); Qs0=Qs0';
			Qsb0=Qsbm(:,:,jj); Qsb0=Qsb0';
			Qsm0=Qsmm(:,:,jj); Qsm0=Qsm0';
			% Convert to GRIDobjs
			QS=GRIDobj(lon(lft_ix:rgt_ix),lat(bot_ix:top_ix),Qs0);
			QSB=GRIDobj(lon(lft_ix:rgt_ix),lat(bot_ix:top_ix),Qsb0);
			QSM=GRIDobj(lon(lft_ix:rgt_ix),lat(bot_ix:top_ix),Qsm0);
			% Resample to DEM
			QSr=resample(QS,DEMLL,'nearest');
			QSBr=resample(QSB,DEMLL,'nearest');
			QSMr=resample(QSM,DEMLL,'nearest');
			QSrr=reproject2utm(QSr,90);
			QSBrr=reproject2utm(QSBr,90);
			QSMrr=reproject2utm(QSMr,90);

			R=QSrr+QSBrr+QSMrr;
			R.Z(R.Z<0)=0;
			RR=flowacc(FD,R)*(DEM.cellsize^2);
			runoff(jj)=(RR.Z(oix)/DA.Z(oix))*(10*100);
		end

		file_out=['routed_gages' filesep 'gage_' shp(ii,1).GAGE_ID]
		T=table(runoff);
		writetable(T,file_out);

	end
end

