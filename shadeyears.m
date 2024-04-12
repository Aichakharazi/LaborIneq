function []=shadeyears(start,finish);

% shadeyears.m  -- feed in your own dates rather than nber
%
%  Routine to shade the dates in a figure

%load nberdatesQ;  % start and finish contain the recession periods, back to 1857
curax=axis;
indx1=find(finish>curax(1));  % First recession to include;
indx2=find(start<curax(2));  % Last recession to include;
indx1=indx1(1);
indx2=indx2(length(indx2));
if start(indx1)<curax(1);
  start(indx1)=curax(1);
end;
if finish(indx2)>curax(2);
  finish(indx2)=curax(2);
end;

colorstr=[242,213,213]/255;

shade(start(indx1:indx2),finish(indx1:indx2),colorstr);
