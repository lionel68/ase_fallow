with mhb_buffer as (
	select routcode, st_buffer(st_centroid(st_transform(geom, 25832)), 1000) as buffergeom
	from public.mhb_plot),
agri as (
select b.routcode, st_union(st_intersection(a.geom, b.buffergeom)) as agrigeom
from dlm1.dlm18_aaa_utm32_gf_v2_s a
join mhb_buffer b
on st_intersects(a.geom, b.buffergeom)
where a.oba_18 in (4101, 4102, 4103, 4104, 4109)
group by b.routcode
order by b.routcode),
agri_diff as (
select a.routcode, st_difference(st_transform(a.agrigeom, 3035), b.uniongeom) as agrigeom
from agri a
join public.mhb_swfforest b
on a.routcode = b.routcode
)
select a.routcode, sum(st_length(st_intersection(a.agrigeom, b.uniongeom))) as edge_m
into public.mhb_edge_atkisswf
from agri_diff a
join public.mhb_swfforest b
on a.routcode = b.routcode
group by a.routcode
order by a.routcode
