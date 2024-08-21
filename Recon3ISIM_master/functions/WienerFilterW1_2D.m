function [wFilter] = WienerFilterW1_2D( sp )
    wFilter.sp=sp;
    wFilter.wDenom=updateCache1(sp);
end