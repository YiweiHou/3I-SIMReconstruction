function [wFilter]=WienerFilterW2_2D( sp )
    wFilter.sp=sp;
    wFilter.wDenom=updateCache(sp);
end