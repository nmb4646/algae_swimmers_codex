function [xd, yd] = repeat_traj_rot4(x, y)
% repeat_traj_rot4_attached
% Attach 4 copies of a trajectory end-to-end, where each copy is rotated
% by 90 degrees relative to the previous one, without flipping.
%
% Inputs:
%   x, y : vectors of equal length
%
% Outputs:
%   xd, yd : concatenated attached trajectory
%
% Notes:
%   - Chirality is preserved because only rotations are used.
%   - The first point of each new segment is placed at the last point
%     of the previous segment.

    x = x(:);
    y = y(:);

    if length(x) ~= length(y)
        error('x and y must have the same length.');
    end

    % Shift original so it starts at the origin
    x0 = x - x(1);
    y0 = y - y(1);

    % Store first rotated segment
    xd = x0;
    yd = y0;

    % Current segment to rotate next
    xc = x0;
    yc = y0;

    for i = 2:4
        % Rotate previous segment by +90 degrees (no reflection)
        xr = -yc;
        yr =  xc;

        % Translate so new segment starts where old one ended
        xr = xr - xr(1) + xd(end);
        yr = yr - yr(1) + yd(end);

        % Append, skipping first point to avoid duplication
        xd = [xd; xr(2:end)];
        yd = [yd; yr(2:end)];

        % Update current segment
        xc = xr - xr(1);
        yc = yr - yr(1);
    end
end