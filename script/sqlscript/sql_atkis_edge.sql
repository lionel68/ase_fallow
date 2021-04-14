create table mhb_1kmbuffer_edge as
select a.routcode, sum(st_length(st_intersection(a.geom, b.agrigeom))) as edge_m
from public.mhb_1kmbuffer_forest a
join public.mhb_1kmbuffer_agri b
on a.routcode = b.routcode
group by a.routcode
order by a.routcode