% clipping and scaling for [-1 1] output range

function [xclipr,xclipi,sc] = clipsc(x,rclip)
xmax = max(max(real(x)),max(imag(x)));
xmin = min(min(real(x)),min(imag(x)));
xamp = xmax - xmin;

xmin = min(xmin + xamp*rclip, -1);
xmax = max(xmax - xamp*rclip, 1);

xclipr = real(x);
xclipr(xclipr<xmin) = xmin;
xclipr(xclipr>xmax) = xmax;

xclipi = imag(x);
xclipi(xclipi<xmin) = xmin;
xclipi(xclipi>xmax) = xmax;

sc = max(abs(xmin),abs(xmax));
xclipr = xclipr/sc;
xclipi = xclipi/sc;
end