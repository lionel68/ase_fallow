select b.routcode, st_union(st_intersection(a.geom, st_buffer(st_transform(st_centroid(b.geom), 25832), 1000))) as agrigeom
from dlm1.dlm18_aaa_utm32_gf_v2_s a
join public.mhb_plot b
on st_intersects(a.geom, st_buffer(st_transform(st_centroid(b.geom), 25832), 1000))
where a.oba_18 in (4101, 4102, 4103, 4104, 4109)
group by b.routcode
order by b.routcode
