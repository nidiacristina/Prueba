function [xo, yo, i_zone] = smpm_utm2ll( xi, yi, i_zone, i_type )
% [xo, yo, i_zone] = smpm_utm2ll( xi, yi, i_zone, i_type )
%
% Converts between (latitude, longitude) to UTM coordinates.
%
% Takes 4 arguments:
%
%  xi     - Vector, of length N, of input coordinates.  Either UTM or
%           latitudes depending on i_type.
%  yi     - Vector, of length N, of input coordinates.  Either UTM or
%           longitudes depending on i_type.
%  i_zone - Specifies the UTM zone to use.  Will be calculated if 0 is
%           provided.  Ignored if i_type does not specify conversion to UTM
%           coordinates.
%  i_type - Integer specifying the type of coordinates in xi and yi.  If
%           specified as 1, xi and yi are latitudes and longitudes,
%           respectively, otherwise xi and yi are UTM coordinates.
%
% Returns 3 values:
%
%  xo     - Vector, of length N, of output coordinates.  Either UTM or
%           latitudes and will be the opposite of xi.
%  yo     - Vector, of length N, of output coordinates.  Either UTM or
%           longitudes and will be the opposite of yi.
%  i_zone - UTM zone used, if UTM coordinates were requested.  If
%           latitude and longitudes were requested, this will remain the
%           same as the input argument.
%

% Originally written in Fortran by Scott Hensley (1993/07/22) and adapted
% to MATLAB by Yuri Fialko (1999/11/25).

% magic constants pass down from upon high.
r_dtor    = 1.74532925199d-2;
i_ft      = 0;
a_griddes = ['C', 'D', 'E', 'F', 'G', 'H', 'J', ...
             'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', ...
             'V', 'W', 'X'];
r_a       = 6378206.4d0;
r_e2      = 0.00676865799761d0;
r_k0      = 0.9996d0;    % scale at center
r_lat0    = 0.0d0;
r_fe      = 5e5;
r_fn(1)   = 0;
r_fn(2)   = 1e7;  % N-S hemispheres??

% begin conversion.
r_ep2  = r_e2 / (1.d0 - r_e2);
r_e4   = r_e2^2;
r_e6   = r_e2^3;
r_dtor = pi/180;

% convert (lat, lon) to UTM.
if i_type == 1

    xi = xi * r_dtor;
    yi = yi * r_dtor;
    if (i_zone == 0)
        i_zone = fix( mod( xi+3.d0*pi, 2.d0*pi ) / (r_dtor * 6.d0) ) + 1;
        i_zone = max( min( i_zone, 60 ), 1 );
    end
    r_lon0 = -pi + 6.d0 * r_dtor * (i_zone - 1) + 3.d0 * r_dtor; % central meridian

    r_n  = r_a / sqrt( 1.d0 - r_e2 * sin( yi )^2 );
    r_t  = tan( yi )^2;
    r_t2 = r_t^2;
    r_c  = r_ep2 * cos( yi )^2;
    r_ba = (xi - r_lon0) * cos( yi );
    r_a2 = r_ba^2;
    r_a3 = r_ba * r_a2;
    r_a4 = r_ba * r_a3;
    r_a5 = r_ba * r_a4;
    r_a6 = r_ba * r_a5;
    r_m  = r_a * ((1.0d0 - r_e2/4 - 3.0d0 * r_e4 / 64.0d0 - 5.d0 * r_e6 / 256.d0) * yi - (3.d0 * r_e2 / 8.d0 + 3.d0 * r_e4 / 32.d0 + 45.d0 * r_e6 / 1024.d0) * sin( 2.d0 * yi ) +  (15.d0 * r_e4 / 256.d0 + 45.d0 * r_e6 / 1024.d0) * sin( 4.d0 * yi ) - (35.d0 * r_e6 / 3072.d0) * sin( 6.d0 * yi ));
    r_m0 = r_a * ((1.d0 - r_e2/4 - 3.d0 * r_e4 / 64.d0 - 5.d0 * r_e6 / 256.d0) * r_lat0 - (3.d0 * r_e2 / 8.d0 + 3.d0 * r_e4 / 32.d0 + 45.d0 * r_e6 / 1024.d0) * sin( 2.d0 * r_lat0 ) + (15.d0 * r_e4 / 256.d0 + 45.d0 * r_e6 / 1024.d0) * sin( 4.d0 * r_lat0 ) - (35.d0 * r_e6 / 3072.d0) * sin( 6.d0 * r_lat0 ));

    r_v(1) = r_k0 * r_n * (r_ba + (1.d0 - r_t + r_c) * r_a3 / 6.d0 + (5.d0 - 18.d0 * r_t + r_t2 + 72.d0 * r_c - 58.d0 * r_ep2) * r_a5 / 120.d0);
    r_v(1) = r_v(1) + r_fe;

    r_v(2) = r_k0 * (r_m - r_m0 + r_n * tan( yi ) * (r_a2 / 2.d0 + (5.d0 - r_t + 9.d0 * r_c + 4.d0 * r_c^2) * (r_a4 / 24.d0) + (61.d0 - 58.d0 * r_t + r_t2 + 600.d0 * r_c - 330.d0 * r_ep2) * (r_a6 / 720.d0) ));
    if yi >= 0
        r_v(2) = r_v(2) + r_fn(1);
    else
        r_v(2) = r_v(2) + r_fn(2);
    end

    r_k = r_k0 * (1.d0 + (1.d0 + r_ep2 * cos( yi )^2) * (r_v(1) - r_fe)^2 / (2.d0 * (r_k0^2) * r_n^2));

    i_gi   = fix( (yi / r_dtor + 80.d0) / 8.d0 ) + 1;
    i_gi   = max( min( i_gi, 20), 1 );
    a_grid = a_griddes(i_gi);

    xo = r_v(1);
    yo = r_v(2);

% convert UTM to (lat, lon).
elseif i_type == 2

    r_vu(1) = xi - r_fe;
    r_vu(2) = yi;
    if r_vu(2) >= r_fn(2)        % LOOK at this if in the S hemisphere!!
        r_vu(2) = yi - r_fn(2);
    end
    r_lon0 = -pi + 6.d0 * r_dtor * (i_zone - 1) + 3.d0 * r_dtor;

    r_et   = sqrt( 1.d0 - r_e2 );
    r_e1   = (1.d0 - r_et) / (1.d0 + r_et);
    r_e12  = r_e1^2;
    r_e13  = r_e1 * r_e12;
    r_e14  = r_e1 * r_e13;
    r_m    = r_vu(2) / r_k0;
    r_mu   = r_m / (r_a * (1.d0 - r_e2 / 4.d0 - 3.d0 * r_e4 / 64.d0 - 5.d0 * r_e6 / 256.d0));
    r_lat1 = r_mu + (3.d0 * r_e1 / 2.d0 - 27.d0 * r_e13 / 32.d0) * sin( 2.d0 * r_mu ) + (21.d0 * r_e12 / 16.d0 - 55.d0 * r_e14 / 32.d0) * sin( 4.d0 * r_mu ) + (51.d0 * r_e13 / 96.d0) * sin( 6.d0 * r_mu ) +(1097.d0 * r_e14 / 512.d0) * sin( 8.d0 * r_mu );

    r_n  = r_a / sqrt( 1.d0 - r_e2 * sin( r_lat1 )^2 );
    r_r  = (r_a * (1.d0 - r_e2)) / sqrt( 1.d0 - r_e2 * sin( r_lat1 )^2 )^3;
    r_t  = tan( r_lat1 )^2;
    r_t2 = r_t^2;
    r_c  = r_ep2 * cos( r_lat1 )^2;
    r_c2 = r_c^2;
    r_d  = r_vu(1) / (r_n * r_k0);
    r_d2 = r_d^2;
    r_d3 = r_d2 * r_d;
    r_d4 = r_d3 * r_d;
    r_d5 = r_d4 * r_d;
    r_d6 = r_d5 * r_d;

    yo = r_lat1 - (r_n * tan( r_lat1 ) / r_r) * (r_d2 / 2.d0 + (5.d0 + 3.d0 * r_t + 10.d0 * r_c - 4.d0 * r_c2 - 9.d0 * r_ep2) * r_d4 / 24.d0 +( 61.d0 + 90 * r_t + 298.d0 * r_c + 45.d0 * r_t2 - 252.d0 * r_ep2 - 3.d0 * r_c2) * (r_d6 / 720.d0));
    xo = r_lon0 + (r_d - (1.d0 + 2.d0 * r_t + r_c) * r_d3 / 6.d0 + (5.d0 - 2.d0 * r_c + 28.d0 * r_t - 3.d0 * r_c2 + 8.d0 * r_ep2 + 24.d0 * r_t2) * (r_d5 / 120.d0)) / cos( r_lat1 );
    xo = xo / r_dtor;
    yo = yo / r_dtor;

end
