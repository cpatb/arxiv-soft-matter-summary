# What is known about schooling like this, and what is not

A reading of the collective-behaviour literature oriented specifically at the question
these four videos can address: a dense, cohesive ball of catfish fry in shallow
freshwater, filmed from above with no calibration.

Everything below that carries a citation was checked against the source during this
session. Where I am relying on background knowledge that I did not verify, I say so
explicitly rather than dressing it up as a citation.

---

## 1. What the animals are

The behaviour in the videos is a well-described life-history stage of North American
bullhead catfish (genus *Ameiurus*). After the yolk sac is absorbed, the fry rise off
the nest and form a dense, dark, ball-shaped school that one or both parents guard and
herd, driving off other fish and physically retrieving strays by mouth. Parental care
lasts up to about a month
([Animal Diversity Web, *Ameiurus nebulosus*](https://animaldiversity.org/accounts/Ameiurus_nebulosus/)).
The Michigan DNR describes the same thing from the shoreline view that matches these
recordings almost exactly: "small clouds of bullhead fry travel the shoreline as the
parents circle the school to keep them together for about 2 weeks"
([Michigan DNR, Bullheads](https://www.michigan.gov/dnr/education/michigan-species/fish-species/bullheads)).

The recordings are from 8 August, in warm shallow nearshore water with abundant weed
growth over a sand and silt bottom. The habitat and season both fit. Two candidate
species, black bullhead (*A. melas*) and brown bullhead (*A. nebulosus*), cannot be
separated from overhead silhouettes at this resolution. **The identification should be
reported as "bullhead catfish fry, *Ameiurus* sp." and not narrowed further.** A single
photograph of an adult, or of a fry in the hand, would settle it.

One caution worth stating plainly. I have not verified from a source that dense dark fry
balls in Lake Michigan nearshore water are always *Ameiurus*. Several other fish produce
dark schooling young. The behavioural match, a compact guarded ball rather than a loose
shoal, is the strongest evidence. It is good evidence, but it is not a specimen.

## 2. The physics tradition: flocks as critical systems

The framing the user has in mind is real and is the right one to reach for.

Giorgio Parisi, with Andrea Cavagna, Irene Giardina and colleagues, measured
three-dimensional positions and velocities of every bird in wild starling flocks over
Rome and asked a statistical-mechanics question of them: how far apart can two birds be
and still have correlated deviations from the flock's mean motion? They found the answer
was *as far as the flock is wide*. The correlation length grew in proportion to flock
size rather than saturating at some fixed interaction range. That is scale-free correlation, a
hallmark of a system sitting near a critical point
([Cavagna et al., PNAS 107:11865, 2010](https://www.pnas.org/content/107/26/11865);
[PubMed 20547832](https://pubmed.ncbi.nlm.nih.gov/20547832/)). The interpretation is that
each bird's *influence* reaches much further than its perception does, which is what lets
a flock respond to a hawk as a single object.

Two companion results matter for how one analyses a school:

- **Interaction is topological, not metric.** Starlings track a roughly fixed *number* of
  neighbours (about seven) rather than everything within a fixed distance, which keeps
  the flock cohesive as density varies (Ballerini et al., PNAS 105:1232, 2008, widely
  cited in the sources above; I did not fetch the original this session).
- **The scale-free property is not automatic.** Cavagna's group has spent a decade
  working out what kind of model reproduces it, converging on "marginal" speed
  confinement, a potential that barely resists small speed deviations but strongly
  suppresses large ones
  ([Cavagna et al. 2021, arXiv:2101.09748](https://arxiv.org/abs/2101.09748);
  [Cavagna et al. 2018, arXiv:1812.07522](https://arxiv.org/abs/1812.07522)).
  Ordinary ferromagnetic alignment models do *not* give scale-free speed correlations in
  their ordered phase. So finding or not finding it in a fish school is informative
  either way.

The same group extended the approach to insects: wild midge swarms show strong
correlation with almost no global order, and finite-size scaling places them near
criticality (Attanasi et al., PRL 113:238102, 2014).

## 3. What has been done with fish

Fish are extremely well studied, but overwhelmingly **in tanks**, and the freshwater
model species is already established.

- **Golden shiners** (*Notemigonus crysoleucas*), a freshwater minnow, are the workhorse
  of the Couzin lab. Katz et al. inferred the effective interaction forces directly from
  trajectories, finding speed regulation rather than explicit alignment to be the dominant
  channel ([PNAS 108:18720, 2011](https://www.pnas.org/content/108/46/18720)). Tunstrøm et
  al. mapped the same species onto order parameters and identified at least three
  dynamically stable collective states: swarm, mill, and polarised
  ([PLoS Comput Biol 9:e1002915, 2013](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1002915)).
  So "freshwater schooling has never been studied" is not correct. It is the standard
  system. What has not been done is the *field* version of it.
- **Criticality in fish** has been probed via turning avalanches, which show critical
  signatures ([Phys Rev Research 6:033270, 2024](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.6.033270)).
- **Environmental modulation**: light level tunes the transition between disordered,
  polarised and milling states ([Lafoux et al., arXiv:2301.09577](https://arxiv.org/abs/2301.09577));
  obstacle density screens social interaction and breaks schooling
  ([Ventéjou et al., arXiv:2402.03123](https://arxiv.org/abs/2402.03123)).
- **Nearest-neighbour distance** is the standard cohesion statistic, alongside speed and
  polarisation; its distribution is typically close to log-normal
  (see the zebrafish length-scale analysis, [arXiv:2405.09469](https://arxiv.org/abs/2405.09469)).
- **Field work exists but is rare.** Ribeiro et al. tracked silver carp in earthen ponds
  and showed shade reduces polarisation and neighbour stability
  ([arXiv:2211.02677](https://arxiv.org/abs/2211.02677)). Freshwater, outdoors and
  quantitative, it is the closest published analogue to what these videos are.

## 4. The methodological gap these videos sit in

Aerial video of wild groups is an active and fast-moving methods area. A 2025 review in
*Methods in Ecology and Evolution* covers drones plus computer vision for collective
animal behaviour
([Kline et al. 2025](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.70128)),
and a 2025 *IJCV* paper builds tracking specifically for drone footage of swarms
([Watching Swarm Dynamics from Above](https://link.springer.com/article/10.1007/s11263-025-02713-y)).

The important admission in that literature, and the one that governs this project: for
dense aggregations the realistic target is often **the group as an object**, its
position, spatial extent and shape, rather than every individual inside it. That is not
a compromise invented here; it is the standard recommendation when silhouettes overlap.

## 5. Where the actual gap is

Putting it together:

| Element | Status |
|---|---|
| Scale-free correlation analysis of a bird flock | Done, foundational (Cavagna/Parisi) |
| Same analysis for midge swarms | Done |
| Freshwater fish schooling, quantitative | Done extensively, but in tanks (golden shiners) |
| Freshwater fish schooling, quantitative, in the field | Rare; silver carp in ponds is the nearest |
| **Ictalurid (catfish) fry schools, quantitative anything** | **I found none** |
| Fission/fusion statistics of a guarded fry ball | I found none |

Juvenile bullhead schooling is documented **descriptively**, the behaviour is in the
natural-history literature and in agency species accounts, traced in these searches to a
1986 reference (Wang) via the Smithsonian NEMESIS
[*Ameiurus melas* summary](https://invasions.si.edu/nemesis/species_summary/164039). What
I could not find, across arXiv and general search, is any quantitative treatment: no
polarisation, no nearest-neighbour distributions, no correlation lengths, no shape
statistics, nothing that treats a bullhead fry ball as a collective-behaviour system.

**A negative literature search is weak evidence.** It shows that the work is not easy to
find, not that it does not exist. Fisheries and ichthyology journals from the 1960s-80s
are poorly indexed and are exactly where descriptive work on bullhead broods would live.
Before any claim of novelty is made in public, someone should search *Copeia*,
*Transactions of the American Fisheries Society*, and the *Journal of Fish Biology*
directly.

With that caveat, the honest positioning is:

> A dense guarded catfish fry ball is a collective-behaviour system that appears not to
> have been measured quantitatively, and it differs from every well-studied system in a
> specific and interesting way: it is a *parentally herded* aggregation. Starling flocks
> and shiner shoals are self-organised among peers. This one has an external agent
> actively pushing strays back in.

That difference is what makes it more than a curiosity. In the Vicsek/Couzin framing,
cohesion emerges from local attraction between equals. Here there is an additional
confining potential imposed from outside by the guarding parent, which is, loosely, an
external field applied to a many-body system. Whether the school's internal statistics
look like a self-organised group or like a confined one is a real question, and it is
askable from overhead video.

## 6. What is realistically measurable from *these* recordings

Grounded in what the data actually support, established in the analysis:

**Solid.** Group-level shape and topology (area, elongation, circularity, boundary
convolution, solidity), all dimensionless or normalised to body length; fission and
fusion events; packing fraction and its profile from rind to core; velocity fields by
optical flow; polarisation; the velocity-fluctuation correlation function and its
correlation length.

**Partial.** Individual positions at the sparse fringe only. Body length from isolated
silhouettes.

**Not possible.** Complete individual tracking. The fish occupy a three-dimensional
volume and the camera sees one projection of it, so silhouettes genuinely superimpose, in the dense core I cannot separate individuals by eye either, at any magnification.
This is an information limit of the recording, not a shortcoming of the detector, and it
should be stated rather than papered over with a tracker that emits confident nonsense.

**Absent entirely.** Metric scale. There is no ruler, no calibration object, and the
camera height is unknown and drifting. Every length in this project is expressed in fish
body lengths, measured from the fish themselves in the same frame, which has the
convenient side effect of cancelling camera-distance drift out of every ratio.

## 7. If this becomes an article

For an urban-ecology or local-natural-history audience, the story that the data can
carry without overreach:

1. A recognisable Chicago place, a specific date, and a behaviour most people walk past.
2. The parent-offspring mechanism, which is genuinely charming and well documented.
3. The physics framing, the same mathematics used on starling murmurations over Rome,
   applied to a catfish ball off Montrose, with the honest note that the fish version is
   two-dimensional and cruder.
4. Numbers that are real: school dimensions in body lengths, how the shape changes,
   how often it splits and rejoins, how far coordination reaches.

What to avoid: any claim that this is the first observation of the behaviour (it is not),
any metric size or speed in centimetres (unmeasurable here), and any individual-level
statistic presented as complete.
