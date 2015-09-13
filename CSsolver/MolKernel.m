function img_kernel = MolKernel(x_inx,y_inx,x_pos,y_pos,mag)
%Simulation
img_kernel = exp(-0.5*(((x_inx-x_pos)/mag).^2+((y_inx-y_pos)/mag).^2))/(sqrt(2*pi)*mag);
           
%110123 A647 data
%img_kernel = 0.7082 * exp(-2*(((x_inx-x_pos)/mag).^2+((y_inx-y_pos)/mag).^2)/1.75345^2) ...
%               +0.3209*exp(-2*(sqrt(((x_inx-x_pos)/mag).^2+((y_inx-y_pos)/mag).^2)-0.43046).^2/1.97672^2);
           
%110211 mEos Data
% img_kernel = 0.87694 * exp(-2*(((x_inx-x_pos)/mag).^2+((y_inx-y_pos)/mag).^2)/1.69217^2) ...
%                + 0.12321 * exp(-2*(((x_inx-x_pos)/mag).^2+((y_inx-y_pos)/mag).^2)/2.38311);
end
