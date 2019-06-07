% clipping and scaling for [-1 1] output range

function [xclipr,xclipi,sc] = clipsc(x,rclip)

xclipr = real(x);
xclipi = imag(x);

xmax = max(max(real(x)),max(imag(x)));
xmin = min(min(real(x)),min(imag(x)));

if rclip ~= 0
    xamp = xmax - xmin;
    xmin = min(xmin + xamp*rclip, -1);
    xmax = max(xmax - xamp*rclip, 1);

    xclipr(xclipr<xmin) = xmin;
    xclipr(xclipr>xmax) = xmax;

    xclipi(xclipi<xmin) = xmin;
    xclipi(xclipi>xmax) = xmax;
end

sc = max(abs(xmin),abs(xmax));
xclipr = xclipr/sc;
xclipi = xclipi/sc;
end