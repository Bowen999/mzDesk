import theme


def render_substructure():
    theme.tab_header(
        "Substructure Search",
        subtitle="Find molecules that contain a given substructure or SMARTS pattern.",
        kicker="search · SMARTS",
    )
    theme.under_construction(
        "Coming soon",
        note="Draw or type a substructure / SMARTS query, screen a compound set for matches "
             "and see the hits highlighted.",
        eta="Planned — next milestone",
    )
