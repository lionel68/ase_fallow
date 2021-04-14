with mhb_buffer as
( select routcode, 
st_buffer(st_centroid(st_transform(geom, 25832)), 2500) as buffergeom
from public.mhb_plot
)
select mhb_buffer.routcode, a.oba_18, b.obb, b.obabez,
sum(st_area(st_intersection(buffergeom, a.geom)))  as atkis_sqm
from dlm1.dlm18_aaa_utm32_gf_v2_s a
join dlm1.atkisgruppierung b
on a.oba_18 = b.oba
join mhb_buffer
on st_intersects(buffergeom, a.geom)
group by mhb_buffer.routcode, a.oba_18, b.obb, b.obabez
order by mhb_buffer.routcode, a.oba_18