from __future__ import annotations

from sirius_prompts import (
    get_core_guardrails,
    get_mode_prompt,
    get_output_contract,
    get_system_sirius_prompt,
)


def build_system_prompt(intent: str, has_docs: bool, has_pubmed: bool) -> str:
    sections = [
        get_core_guardrails(),
        get_system_sirius_prompt(),
        get_mode_prompt(intent),
        get_output_contract(intent),
        (
            "RUNTIME CONTEXT FLAGS\n"
            f"- has_workdir_documents: {'yes' if has_docs else 'no'}\n"
            f"- has_pubmed_context: {'yes' if has_pubmed else 'no'}"
        ),
    ]
    return "\n\n".join([s for s in sections if s])
