function [subArrival] = subset_tlArrival_tlPicker(tlArrival, aPhase, aStation, aEvent, aEventType)
%SUBSET_tlARRIVAL -- Choose a subset of tlArrival
%                    (TomoLab Utility)
%
%   [subArrival] = subset_tlArrival_tlPicker(tlArrival, aPhase, aStation, aEvent, aEventType)
%   chooses a subset of the tlArrival structure.  Subset is chosen given
%   one or all of the calling arguments aPhase, aStation, aEvent or aEventType.
%   Arguments left empty are not used in making the subset.
%   Multiple phases or stations are entered as cell arrays
%
%
%   Returns empty fields if no matches.
%   Otherwise returns:
%
%          subArrival.eventid       see tlArrival
%          subArrival.phase         see tlArrival
%          subArrival.error         see tlArrival
%          subArrival.station       see tlArrival
%          subArrival.time          see tlArrival
%          subArrival.eventtype     see tlArrival
%          subArrival.num           number of arrivals returned
%          subArrival.ptr           pointer to indices of tlArrival that
%                                   comprise subset
%
% Modified by William Wilcock to accept vector inputs and make use of
% functions cell2string and vector_match

%  Copyright 2010 Blue Tech Seismics, Inc.

subArrival.eventid   = [];
subArrival.phase     = [];
subArrival.error     = [];
subArrival.station   = [];
subArrival.time      = [];
subArrival.num       = [];
subArrival.eventtype      = [];
subArrival.phase     = {};
subArrival.nphases   = [];

if nargin ~= 5
    error('5 input arguments: subset_tlArrival_tlPicker(tlArrival, aPhase, aStation, aEvent, aEventType)')
end

if isempty(aPhase) && isempty(aStation) && isempty(aEvent) && isempty(aEventType)
    error('subset_tlArrival_tlPicker:  One of the last 4 input arguments must be non-empty')
end

% Make phase and station strings single element cell structures 
if ischar(aPhase)
  aPhase = {aPhase};
end
if ischar(aStation)
  aStation = {aStation};
end

J=1:length(tlArrival.eventid);

% Subset for station(s)

if ~isempty(aStation)
    if ~all(cellfun('isclass',aStation,'char'))
        error('aStation must be character string or cell structure of character strings')
    end
    jsta = false(size(tlArrival.station(J)));
    for i = 1:length(aStation)
      jsta = jsta | strcmp(aStation(i), tlArrival.station(J));
    end
    if any(jsta)
        J = J(jsta);
    else
        warning(['Skipping station(s) (', cell2string(aPhase), ') No arrivals (subset_tlArrival_tlPicker)'])
        return
    end
end

% Subset for event(s)

if ~isempty(aEvent)
    if ~isnumeric(aEvent)
        error('aEvent must be numeric')
    end
    jevt = vector_match(tlArrival.eventid(J),aEvent);
    if any(jevt)
        J = J(jevt);
    else
        warning(['Skipping event(s) (', int2str(aEvent),') No arrivals (subset_tlArrival_tlPicker)'])
        return
    end
end

%  Subset for phase(s)

if ~isempty(aPhase)
    if ~all(cellfun('isclass',aPhase,'char'))
        error('aPhase must be character string or cell structure of character strings')
    end
    jphs = false(size(tlArrival.phase(J)));
    for i=1:length(aPhase)
      jphs = jphs | strcmp(aPhase(i), tlArrival.phase(J));
    end
    if any(jphs)
        J = J(jphs);
    else
        warning(['Skipping phase(s) (', cell2string(aPhase),') No arrivals (subset_tlArrival_tlPicker)'])
        return
    end
end

%  Subset for event type

if ~isempty(aEventType)
    if ~isnumeric(aEventType)
        error('aEventType must be numeric')
    end
    jtyp = vector_match(tlArrival.eventtype(J),aEventType);
    if any(jtyp)
        J = J(jtyp);
    else
        warning(['Skipping eventtype(s) (', int2str(aEventType),') No arrivals (subset_tlArrival_tlPicker)'])
        return
    end
end

%  Station-phase pair has data, so set values.

subArrival.eventid   = tlArrival.eventid(J);
subArrival.phase     = tlArrival.phase(J);
subArrival.error     = tlArrival.error(J);
subArrival.station   = tlArrival.station(J);
subArrival.time      = tlArrival.time(J);
subArrival.eventtype = tlArrival.eventtype(J);
subArrival.num       = length(subArrival.time);
subArrival.ptr       = J(:);
subArrival.phaselist = unique(tlArrival.phase(J));
subArrival.nphases   = length(subArrival.phaselist);


