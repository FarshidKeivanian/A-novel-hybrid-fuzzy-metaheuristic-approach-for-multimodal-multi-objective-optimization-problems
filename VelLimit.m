function [VelMin, VelMax] =VelLimit(GlobalBestPosition, xPosition, it, up, down, Alpha)

VelMax = +Alpha.*(((up-down)./up).*abs((GlobalBestPosition - xPosition)/it));
VelMin = -(VelMax);

end