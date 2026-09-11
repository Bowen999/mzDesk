import streamlit as st
import theme


def render_substructure():
    theme.tab_header(
        "Substructure Search",
        subtitle="Find molecules that contain a given substructure or SMARTS pattern.",
        kicker="search · SMARTS",
    )
    theme.under_construction(
        "Substructure Search",
        note="This module will let you draw or type a substructure / SMARTS query and screen "
             "a compound set for matches, with highlighted hits.",
        eta="Planned — next milestone",
    )
