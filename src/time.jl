
function porewater(p::Proposal, sc::SedimentColumn, k::NamedTuple, climhist::NamedTuple, seawater::Seawater, ka_dt::Int)

    isd = searchsortedfirst(climhist.t, p.onset, rev=true)

    isd = ifelse(isd<climhist.n, isd, climhist.n)

    @inbounds for t = isd:climhist.n

        PorewaterDiffusion.boundaryconditions(sc.Cl.o[1], sc.O.o[1], climhist.x[t], p.sea2frz, p.frz2mlt, p.dmlt, p.dfrz, seawater.Cl, seawater.O, k.dz, k.dt)

        @inbounds for j = Base.OneTo(ka_dt)
            diffuseadvectcolumn!(sc,k)
        end
    end
end