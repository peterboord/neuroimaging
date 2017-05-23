
function [data, index]=convertTDMS(filename);

%convertTDMS - function to convert Labview TDMS data files into .mat files.
%   If called with no input, user selects from file open dialog box.  A
%   .mat file with the same filename as the TDMS file is automatically
%   written to the same directory (warning - will over write .mat file of
%   the same name.
%
%   TDMS format is based on information provided by National Instruments
%   at:    http://zone.ni.com/devzone/cda/tut/p/id/5696
%
% convertTDMS(filename)
%
%       Inputs:
%               filename - filename to be converted.  If not supplied, user
%                 is provided dialog box to open file.  Can be a cell array
%                 of files for bulk conversion.
%
%       Outputs:
%               ob - structure with all of the data objects
%               index - index of the information in ob
%

%---------------------------------------------
%Brad Humphreys - v1.0 2008-04-23
%ZIN Technologies
%---------------------------------------------
%---------------------------------------------
%Brad Humphreys - v1.1 2008-07-3
%ZIN Technologies
%-Added abilty for timestamp to be a raw data type, not just meta data.
%-Addressed an issue with having a default nsmaples entry for new objects.
%-Added Error trap if file name not found.
%-Corrected significant problem where it was assumed that once an object
%    existsed, it would in in every subsequent segement.  This is not true.
%---------------------------------------------

%---------------------------------------------
%Grant Lohsen - v1.2 2009-11-15
%Georgia Tech Research Institute
%-Converts TDMS v2 files
%Folks, it's not pretty but I don't have time to make it pretty. Enjoy.
%---------------------------------------------


%---------------------------------------------
%Jeff Sitterle - v1.3 2010-01-10
%Georgia Tech Research Institute
%Modified to return all information stored in the TDMS file to inlcude
%name, start time, start time offset, samples per read, total samples, unit
%description, and unit string.  Also provides event time and event
%description in text form
%Vast speed improvement as save was the previous longest task
%---------------------------------------------

%---------------------------------------------
%Grant Lohsen - v1.4 2009-04-15
%Georgia Tech Research Institute
%Reads file header info and stores in the Root Structure.
%---------------------------------------------


dispPerComp = 0;   %set to 1 to display percentage complete


conver='1.4';    %conversion version number

%startingdir=cd; %Get the starting directory

if nargin==0                %If file is not provided,prompt user
    [filename,pathname,filterindex] = uigetfile({'*.*',  'All Files (*.*)'},'Choose File');
    filename=fullfile(pathname,filename);
end

%Create filelist
if iscell(filename)         %If it is a group of files
    numfiles=max(size(filename));
    infilename=filename;
else
    numfiles=1;             %If only one file name is provided
    infilename{1}=filename;
end



for fnum=1:numfiles
    disp([datestr(now,13) ' Begining conversion of:  '  infilename{fnum}])
    
    bytesInfile = dir(infilename{fnum});    %Get size of file.  Needed for later estimation of variable size.
    
    if isempty(bytesInfile)
        error(['Error: ' infilename{fnum} ' not Found']);
    end
    filesize = bytesInfile.bytes;
    
    %Initialize variables for each file conversion
    index.names=[];
    ob=[];
    %firstRawData=1;
    segCnt=0;
    
    fid=fopen(infilename{fnum});
    rawDataInThisSeg = [];
    fseek(fid, 0, 1);
    eoff = ftell(fid);
    fseek(fid, 0, -1);
    perComp = 0;
    while (ftell(fid) ~= eoff)
 %   while ~feof(fid) %does not work
        if dispPerComp == 1
            if (100*ftell(fid)/eoff) > (perComp + 10) 
                disp(sprintf('Percent Complete: %2.1f%%', (100*ftell(fid)/eoff)))
                perComp = (100*ftell(fid)/eoff);
            end
        end
        Ttag=fread(fid,1,'uint8');
        Dtag=fread(fid,1,'uint8');
        Stag=fread(fid,1,'uint8');
        mtag=fread(fid,1,'uint8');
        if Ttag==84 & Dtag==68 & Stag==83 & mtag==109
            segCnt=segCnt+1;
            if ~isempty(rawDataInThisSeg)
            rawDataInThisSeg(1:length(rawDataInThisSeg))=0;      % Reset the raw data indicator
            end
            
            %ToC Field
            Toc=fread(fid,1,'uint32');
            kTocMetaData=bitget(Toc,2);
            kTocNewObject=bitget(Toc,3);
            kTocRawData=bitget(Toc,4);
            
            %Segment
            vernum=fread(fid,1,'uint32');
            if vernum ~= 4713 
                    sprintf('Warning: This code is designed for version 4713 of the labview TDMS format.\nYou are attempting to import a file from version %s . Aborting Now', num2str(vernum))
                return
            end
            
            segLength=fread(fid,1,'uint64');
            
            metaRawLength=fread(fid,1,'uint64');
            
            %% Process Meta Data
            if kTocMetaData                                     %If there is meta data in this segment
                numNewObjInSeg=fread(fid,1,'uint32');
                
                for q=1:numNewObjInSeg
                    
                    obLength=fread(fid,1,'uint32');             %Get the length of the objects name
                    obname=[fread(fid,obLength,'uint8=>char')]';%Get the objects name
                    
                    %Fix Object Name
                    if strcmp(obname,'/')
                        obname='Root';
                    else
                        obname=fixcharformatlab(obname);
                    end
                    
                    if ~isfield(ob,obname)                         %If the object does not already exist, create it
                        index.names{end+1}=obname;
                        ob.(obname)=[];
                        obnum=max(size(index.names));               %Get the object number
                        newob(obnum)=1;                             %Brand new object
                    else                                            %if it does exist, get it's index and object number
                        newob(obnum)=0;
                        obnum=find(strcmp(index.names,obname)==1,1,'last');
                    end
                    
                    rawdataindex=fread(fid,1,'uint32');             %Get the raw data Index
                    if rawdataindex==0                              %No raw data assigned to this object in this segment
                        index.entry(obnum)=0;
                        index.dataType(obnum)=1;
                        index.arrayDim(obnum)=0;
                        index.nValues(obnum)=0;
                        index.byteSize(obnum)=0;
                        rawDataInThisSeg(obnum)=0;
                    elseif rawdataindex+1==2^32                     %Objects raw data index matches previous index - no changes
                        if ~isfield(index,'entry')                   %The root object will always have a FFFFFFFF entry
                            rawDataInThisSeg=0;
                        else
                            if max(size(index.entry))<obnum         %In case an object besides the root is added that has no dat but
                                rawDataInThisSeg=0;                    %reports using previous
                            end
                        end
                    else                                            %Get new object information
                        index.entry(obnum)=rawdataindex;
                        index.dataType(obnum)=fread(fid,1,'uint32');
                        index.arrayDim(obnum)=fread(fid,1,'uint32');
%                         disp('nvals');
%                         tango = ftell(fid)
%                         if tango == 1633574
%                         disp('halto')
%                         end
                        
                        index.nValues(obnum)=fread(fid,1,'uint64');
                        if index.dataType(obnum)==32                %If the datatype is a string
                            index.byteSize(obnum)=fread(fid,1,'uint64');
                        else
                            index.byteSize(obnum)=0;
                        end
                        rawDataInThisSeg(obnum)=1;
                    end
                    
                    
                    numProps=fread(fid,1,'uint32');
                    for p=1:numProps
                        propNameLength=fread(fid,1,'uint32');
                        propsName=[fread(fid,propNameLength,'uint8=>char')]';
                        propsName=fixcharformatlab(propsName);
                        propsDataType=fread(fid,1,'uint32');
                        propExists=isfield(ob.(obname),propsName);
                        dataExists=isfield(ob.(obname),'data');
                        
                        if dataExists                                               %Get number of data samples for the object in this segment
                            %nsamps=max(size(ob.(obname).data));
                            
                            
                            nsamps=ob.(cname).nsamples+1;
                            
                        else
                            nsamps=0;
                            estNumSeg=floor(filesize/(20+segLength)*1.2);            %Estimate # of Segements.  20 is the number of bytes prior to the segLength Read
                        end
                        
                        if propsDataType==32                                         %If the data type is a string
                            propsValueLength=fread(fid,1,'uint32');
                            propsValue=fread(fid,propsValueLength,'uint8=>char');
                            if propExists
                                cnt=ob.(obname).(propsName).cnt+1;
                                ob.(obname).(propsName).cnt=cnt;
                                ob.(obname).(propsName).value{cnt}=propsValue;
                                ob.(obname).(propsName).samples(cnt)=nsamps;
                            else
                                if strcmp(obname, 'Root')
                                    %header data
                                    %ob.(obname).(propsName).cnt=1;
                                    
                                    ob.(obname).(propsName) = propsValue'; %fixes the annoying habit of the char array being vertical 
                                    %ob.(obname).(propsName).samples(1)=nsamps;

                                else
                                    %group.channel data
                                    ob.(obname).(propsName).cnt=1;
                                    ob.(obname).(propsName).value{estNumSeg}=NaN;
                                    ob.(obname).(propsName).samples(estNumSeg)=0;
                                    ob.(obname).(propsName).value{1}=propsValue;
                                    ob.(obname).(propsName).samples(1)=nsamps;
                                end
                            end
                        else                                                        %Numeric Data type
                            if propsDataType==68                                     %If the data type is a timestamp
                                tsec=fread(fid,1,'uint64')/2^64+fread(fid,1,'uint64');   %time since Jan-1-1904 in seconds
                                propsValue=tsec/86400+695422-5/24;                   %/864000 convert to days; +695422 days from Jan-0-0000 to Jan-1-1904
                            else
                                matType=LV2MatlabDataType(propsDataType);
                                propsValue=fread(fid,1,matType);
                            end
                            if propExists
                                cnt=ob.(obname).(propsName).cnt+1;
                                ob.(obname).(propsName).cnt=cnt;
                                ob.(obname).(propsName).value(cnt)=propsValue;
                                ob.(obname).(propsName).samples(cnt)=nsamps;
                            else
                                ob.(obname).(propsName).cnt=1;
                                ob.(obname).(propsName).value(estNumSeg)=NaN;
                                ob.(obname).(propsName).samples(estNumSeg)=0;
                                ob.(obname).(propsName).value(1)=propsValue;
                                ob.(obname).(propsName).samples(1)=nsamps;
                            end
                        end
                        
                    end
                    
                end
            end
            
            %% Process Raw Data
            
            
            
            
        else
            if ~(isempty(Ttag) & isempty(Dtag) & isempty(Stag) & isempty(mtag))  %On last read, all will be empty
                %error('Unexpected bit stream in file: %s  at segment %d (dec: %d %d %d %d)',infilename{fnum}, segCnt+1, Ttag, Dtag, Stag, mtag);
                fseek(fid, -4, 0);
                %disp('set');
                %ftell(fid)
            end
        end
        
        
        
        for r=1:max(size(index.names))                                      %Loop through the index
            
            if newob(r) & rawDataInThisSeg(r)                               %If brand new object and new raw data, preallocate data arrays
                %                     firstRawData=0;
                newob(r)=0;
                estNumSeg=filesize/(20+segLength);                              %Estimate # of Segements.  20 is the number of bytes prior to the segLength Read
                
                %for b=1:max(size(index.names))
                cname=cell2mat(index.names(r));
                nEstSamples=floor(1.2*estNumSeg*index.nValues(r));
                if index.dataType(r)~=32 | index.dataType(r)~=68             %If the data is numeric type
                    ob.(cname).data(1:(nEstSamples),1)=NaN;
                else                                                         %If the data is string type
                    ob.(cname).data{nEstSamples}=NaN;
                end
                ob.(cname).nsamples=0;
                %end
            end
            if rawDataInThisSeg(r)
                nvals=index.nValues(r);
                if index.dataType(r)~=68                                    %If not a timestamp data type
                    [data cnt]=fread(fid,nvals,LV2MatlabDataType(index.dataType(r)));
                else                                                        %If a timestamp
                    clear data
                    for dcnt=1:nvals
                        tsec=fread(fid,1,'uint64')/2^64+fread(fid,1,'uint64');   %time since Jan-1-1904 in seconds
                        data(1,dcnt)=tsec/86400+695422-5/24;
                    end
                    cnt=dcnt;
                end
                cname=cell2mat(index.names(r));
                if isfield(ob.(cname),'nsamples')
                    ssamples=ob.(cname).nsamples;
                else
                    ssamples=0;
                end
                
                %                     ssamples=ob.(cname).nsamples;
                if index.dataType(r)~=32                                    %If the data is numeric type
                        ob.(cname).data(ssamples+1:ssamples+cnt,1)=data;    
                else                                                        %If the data is string type
                    
                    ob.(cname).data(ssamples+1:ssamples+cnt,1)=data;
                    
                end
                ob.(cname).nsamples=ssamples+cnt;
            end
        end
        %rawDataInThisSeg(1:end)=0;      % Reset the raw data indicator
    end
end


fclose(fid);

%% Clean up preallocated arrays   (preallocation required for speed)
for y=1:max(size(index.names))
    cname=cell2mat(index.names(y));
    if isfield(ob.(cname),'nsamples')
        nsamples=ob.(cname).nsamples;
        if nsamples>0       %Remove any excess from preallocation of data
            if index.dataType(y)~=32
                ob.(cname).data=ob.(cname).data(1:nsamples,1);                  %If the data is numeric type
            else
                ob.(cname).data=ob.(cname).data(1:nsamples,1);                  %If the data is string type
            end
            
            proplist=fieldnames(ob.(cname));    %Remove any excess from preallocation of properties
            for isaac=1:size(proplist,1)
                if isfield(ob.(cname).(proplist{isaac}),'cnt')
                    cnt=ob.(cname).(proplist{isaac}).cnt;
                    ob.(cname).(proplist{isaac}).value=ob.(cname).(proplist{isaac}).value(1:cnt);
                    ob.(cname).(proplist{isaac}).samples=ob.(cname).(proplist{isaac}).samples(1:cnt);
                    ob.(cname).(proplist{isaac})=rmfield(ob.(cname).(proplist{isaac}),'cnt');
                end
            end
            
        end
    end
end



ob.index=index;
ob.conver=conver;
data = postProcess(ob);
% cd(startingdir);
end

function DataStructure = postProcess(ob)

%Modified to return all information stored in the TDMS file to inlcude
%name, start time, start time offset, samples per read, total samples, unit
%description, and unit string.  Also provides event time and event
%description in text form

%Jeff Sitterle
%January 10, 2010

DataStructure.Root = [];
DataStructure.MeasuredData.Name = [];
DataStructure.MeasuredData.Data = [];
DataStructure.Events.Name = [];
DataStructure.Events.Data = [];
%DataStructure.RawOutput = ob;
varNameMask = '';
cntData = 1;
cntEvent = 1;

for i = 1:length(ob.index.names)
    cname=cell2mat(ob.index.names(i));
    if strcmp(cname, 'Root')
        DataStructure.Root = ob.(cname);
        
    end
    if isfield(ob.(cname),'data')
        if strcmp(varNameMask, 'Events')
            DataStructure.Events(cntEvent).Name = char(regexprep(ob.index.names(i), varNameMask, '', 'ignorecase'));
            
            if strcmp(DataStructure.Events(cntEvent).Name, 'Description')
                event_string = char(ob.(cname).data');
                seperator = event_string(1:4);
                locations = findstr(seperator, event_string);
                num_events = max(size(locations));
                for j = 1:num_events
                    if j < num_events
                        DataStructure.Events(cntEvent).Data(j,:) = cellstr(event_string(locations(j)+4:locations(j+1)-1));
                    else
                        DataStructure.Events(cntEvent).Data(j,:) = cellstr(event_string(locations(j)+4:max(size(event_string))));
                    end
                end
            else
                DataStructure.Events(cntEvent).Data = ob.(cname).data;
            end
            cntEvent = cntEvent + 1;
        
        else
            
            DataStructure.MeasuredData(cntData).Name = char(regexprep(ob.index.names(i), varNameMask, '', 'ignorecase'));
            DataStructure.MeasuredData(cntData).Start_Time = ob.(cname).wf_start_time.value;
            DataStructure.MeasuredData(cntData).Start_Time_Offset = ob.(cname).wf_start_offset.value;
            DataStructure.MeasuredData(cntData).Sample_Rate = ob.(cname).wf_increment.value;
            DataStructure.MeasuredData(cntData).Samples_Per_Read = ob.(cname).wf_samples.value;
            DataStructure.MeasuredData(cntData).Total_Samples = ob.(cname).nsamples;
            DataStructure.MeasuredData(cntData).Units_Decription = char(ob.(cname).NI_UnitDescription.value)';
            DataStructure.MeasuredData(cntData).Unit_String = char(ob.(cname).unit_string.value)';
            DataStructure.MeasuredData(cntData).Data = ob.(cname).data;
            cntData = cntData + 1;
        
        end
        
        
    else
        varNameMask = ob.index.names(i);
        varNameMask = varNameMask{:};
    end
    
end


end


function  fixedtext=fixcharformatlab(textin)
%Private Function to remove all text that is not MATLAB variable name
%compatible
textin=strrep(textin,'''','');
textin=strrep(textin,'\','');
textin=strrep(textin,'/Untitled/','');
textin=strrep(textin,'/','.');
textin=strrep(textin,'-','');
textin=strrep(textin,'?','');
textin=strrep(textin,' ','_');
textin=strrep(textin,'.','');
textin=strrep(textin,'[','_');
textin=strrep(textin,']','');
textin=strrep(textin,'%','');
textin=strrep(textin,'#','');
textin=strrep(textin,'(','');
textin=strrep(textin,')','');
textin=strrep(textin,':','');
fixedtext=textin;

end

function matType=LV2MatlabDataType(LVType)
%Cross Refernce Labview TDMS Data type to MATLAB

switch LVType
    case 1   %tdsTypeVoid
        matType='';
    case 2   %tdsTypeI8
        matType='int8';
    case 3   %tdsTypeI16
        matType='int32';
    case 4   %tdsTypeI32
        matType='int32';
    case 5   %tdsTypeI64
        matType='int64';
    case 6   %tdsTypeU8
        matType='uint8';
    case 7   %tdsTypeU16
        matType='uint16';
    case 8   %tdsTypeU32
        matType='uint32';
    case 9   %tdsTypeU64
        matType='uint64';
    case 10  %tdsTypeSingleFloat
        matType='float64';
    case 11  %tdsTypeDoubleFloat
        matType='float64';
    case 12  %tdsTypeExtendedFloat
        matType='';
    case 32  %tdsTypeString
        matType='char';
    case 33  %tdsTypeBoolean
        matType='bit1';
    case 68  %tdsTypeTimeStamp
        matType='bit224';
    otherwise
        matType='';
end

end

