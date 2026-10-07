"""SNP allele and local phaselet anchor used by contact extraction."""
from dataclasses import dataclass


@dataclass(frozen=True)
class ComponentSite:
    ref: str
    alt: str
    delta: int
    component: str
    snp_id: str
