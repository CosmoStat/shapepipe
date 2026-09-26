"""Fingerprint rendered prologues and shells without checkout or tmp paths."""

import hashlib
import json

from tests.workflow.harness import REPO


def params_pin(dag):
    """Return SHA-256 digests for unit_pre and every declared rule's shells.

    Each rule includes all resolved wildcard instances, its shell template,
    and its rendered ``params.pre`` (null for rules without that parameter).
    The normalized payload is saved beside the fixture for review on failure.
    Script hashes and other non-pre params are outside this pin's scope.
    """
    campaign = dag.campaign

    def normalize(value):
        if value is None:
            return None
        return (value.replace(str(campaign.root), "<CAMPAIGN_ROOT>")
                .replace(str(REPO), "<REPO>"))

    def digest(value):
        text = json.dumps(value, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(text.encode()).hexdigest()

    namespace = dag.namespace
    unit_pre = {}
    for stage, (level, _) in sorted(namespace["STAGE_DIR"].items()):
        unit = next(iter(campaign.ready)) if level == "tile" else "2243881"
        unit_pre[stage] = normalize(namespace["unit_pre"](stage, unit))
    rules = {}
    for rule in sorted(dag.workflow.rules, key=lambda rule: rule.name):
        jobs = dag.jobs_for(rule.name)
        # Aggregation-only targets have no shell or pre to expand.
        assert jobs or (not rule.shellcmd and not rule.params), rule.name
        rules[rule.name] = {
            "shell_template": normalize(rule.shellcmd),
            "jobs": [{
                "wildcards": dict(sorted(job.wildcards_dict.items())),
                "pre": normalize(getattr(job.params, "pre", None)),
                "shell": normalize(job.shellcmd),
            } for job in jobs],
        }
    payload = {"unit_pre": unit_pre, "rules": rules}
    (campaign.root / "params_rendered.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )
    return {
        "schema": 1,
        "algorithm": "sha256",
        "sha256": digest(payload),
        "unit_pre": {stage: digest(pre) for stage, pre in unit_pre.items()},
        "rules": {name: digest(rule) for name, rule in rules.items()},
    }
