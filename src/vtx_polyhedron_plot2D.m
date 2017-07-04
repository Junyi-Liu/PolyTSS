function [  ] = vtx_polyhedron_plot2D( vtx, s, n, s_val )
%vtx_polyhedron_plot Summary of this function goes here
%   plot the polyhedron constructed by the input vertices
%   NOTE: input vtx is column-wise

% substitute all symbolic variables by s_val
% vtx_uni is row-wise
vtx_uni = double(subs(vtx', s, s_val*ones(n,1)));

P = Polyhedron('V', vtx_uni);

P.plot;

end

