% clipping for quad MZM, thresh is [min max] for output range

function Eclip = clip(Ein,thresh)
rclip = real(Ein);
rclip(rclip > thresh(2)) = thresh(2);
rclip(rclip < thresh(1)) = thresh(1);

iclip = imag(Ein);
iclip(iclip > thresh(2)) = thresh(2);
iclip(iclip < thresh(1)) = thresh(1);

Eclip = rclip + 1j*iclip;
