insert into public.mhb_bkr
(with bkr2 as
(select bkr_2stell, st_transform(st_union(geom), 4326) geom
from public.boden_klima_raueme
group by bkr_2stell)
select a.routcode, bkr2.bkr_2stell
from public.mhb_plot a
join bkr2
on st_intersects(st_centroid(a.geom), bkr2.geom)
order by a.routcode)