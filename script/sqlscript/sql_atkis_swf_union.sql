with mhb_buffer as (
	select routcode, st_buffer(st_centroid(st_transform(geom, 25832)), 1000) as buffergeom25832,
	st_buffer(st_centroid(st_transform(geom, 3035)), 1000) as buffergeom3035
	from public.mhb_plot
),
forest as (
select b.routcode, st_union(st_intersection(a.geom, b.buffergeom25832)) as forestgeom
from dlm1.dlm18_aaa_utm32_gf_v2_s a
join mhb_buffer b
on st_intersects(a.geom, b.buffergeom25832)
where a.oba_18 in (4107, 4108)
group by b.routcode
),
swf as (
select b.routcode, st_union(st_intersection(a.geom, b.buffergeom3035)) as swfgeom
from copernicus_hrl.swf_2015_vec_de_3035_full a
join mhb_buffer b
on st_intersects(a.geom, b.buffergeom3035)
group by b.routcode
)
select coalesce(a.routcode, b.routcode) as routcode,
st_union(a.swfgeom, st_transform(b.forestgeom, 3035)) as uniongeom
into public.mhb_swfforest
from swf a
full join forest b
on a.routcode = b.routcode
order by routcode