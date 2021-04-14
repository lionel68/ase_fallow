with muni_agri as (select b.reg, st_union(st_intersection(a.geom, st_transform(b.geom, 25832))) as agrigeom
from dlm1.dlm18_aaa_utm32_gf_v2_s a
join public.municipality b
on st_intersects(a.geom, st_transform(b.geom, 25832))
where a.oba_18 in (4101, 4102, 4103, 4104, 4109)
group by b.reg
order by b.reg),
mhb_buffer as (
	select routcode, st_buffer(st_centroid(st_transform(geom, 25832)), 1000) as buffergeom
	from public.mhb_plot),
mhb_muni_agri as (
select b.routcode, a.reg, st_area(st_intersection(a.agrigeom, b.buffergeom)) / st_area(a.agrigeom) as agriprop
from muni_agri a
join mhb_buffer b
on st_intersects(a.agrigeom, b.buffergeom)
)
select b.routcode, a.year, sum(b.agriprop * cast(a.value as decimal(10, 3))) as fallow_ha
into public.mhb_fallow_atkis
from public.agraratlas_2003_2016 a
join mhb_muni_agri b
on a.reg = b.reg
where a.variable in ('FALL', 'SETA')
group by b.routcode, a.year
order by b.routcode, a.year

